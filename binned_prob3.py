#/usr/bin/python
import sys
import numpy as np
from scipy.interpolate import RectBivariateSpline


##########################################################################
# Elim's prob3 interface
##########################################################################
from icecube import prob3 as BargerPropagator
earth_model = "PREM_60layer.dat"
detector_depth = 2.
prop_height   = 20.
barger_prop = BargerPropagator.BargerPropagator(earth_model, detector_depth)
barger_prop.UseMassEigenstates(False)
barger_prop.SetOneMassScaleMode(False)
barger_prop.SetWarningSuppression(True)

def getProb3 (e, cz, pdg, params):
    NuE = 1 ; NuMu = 2 ; NuTau = 3
    sin_sq_theta12 = np.sin(params['theta12'])**2
    sin_sq_theta23 = np.sin(params['theta23'])**2
    sin_sq_theta13 = np.sin(params['theta13'])**2
    dm32 = params['dm31'] - params['dm21']
    kSquared = True
    kNuType = np.array(np.sign(pdg), dtype=np.int)
    out_nu = (  NuE*(np.abs(pdg) == 12) + NuMu*(np.abs(pdg) == 14) + NuTau*(np.abs(pdg) == 16) )
    osc_prob = np.zeros([len(e), 2])
    for i in range(len(e)):
        barger_prop.SetMNS(sin_sq_theta12, sin_sq_theta13, sin_sq_theta23,
                           params['dm21'], dm32, params['deltacp'],
                           float(e[i]), kSquared, int(kNuType[i]))
        barger_prop.DefinePath(float(cz[i]), prop_height)
        barger_prop.propagate(int(kNuType[i]))
        osc_prob[i,:] = [barger_prop.GetProb(NuE, int(out_nu[i])),
                         barger_prop.GetProb(NuMu, int(out_nu[i]))]
    return osc_prob



##########################################################################
# A binned prob3 calculator
##########################################################################
class binned_prob3:
    ebins, czbins, nubar = None, None, None
    emap, czmap, nu_nubar_map = None, None, None
    pmap_from_e, pmap_from_mu = None, None

    pspline = { -12:{}, -14:{}, 
                 12:{}, 14:{} }
    
    oscil_params = {"dm21": 7.49e-5,
                    "dm31": 2.526e-3,
                    "theta23": np.arcsin( np.sqrt(0.440) ),
                    "theta12": np.arcsin( np.sqrt(0.308) ),
                    "theta13": np.arcsin( np.sqrt(0.02163) ),
                    "deltacp": 0,
                    }    
    changed = True
    

    #-----------------------------------------------------------------------
    # Set this sucker up and get the maps that we'll need
    #-----------------------------------------------------------------------
    def __init__(self,
                 ebins = np.logspace(0,3,150), 
                 czbins = np.linspace(-1,1,150)):
        self.ebins = ebins
        self.czbins = czbins
        self.nubar = np.array([-1,1])

        # Make sure to include a nu/nubar part!
        self.emap, self.czmap, self.nu_nubar_map = np.meshgrid(self.ebins, self.czbins, self.nubar, 
                                                               indexing='ij')

        self.pmap_from_e = {12: np.zeros_like(self.emap), 
                            14: np.zeros_like(self.emap), 
                            16: np.zeros_like(self.emap)}

        self.pmap_from_mu = {12: np.zeros_like(self.emap), 
                             14: np.zeros_like(self.emap), 
                             16: np.zeros_like(self.emap)}
    
        self.setParameters()


    #-----------------------------------------------------------------------
    # Set the current oscillation parameters. If these have changed, note
    # that so that we can calculate new oscillation probabilities
    #-----------------------------------------------------------------------
    def setParameters(self,
                      dm21 = 7.49e-5,
                      dm31 = 2.526e-3,
                      theta23 = np.arcsin( np.sqrt(0.440) ),
                      theta12 = np.arcsin( np.sqrt(0.308) ),
                      theta13 = np.arcsin( np.sqrt(0.02163) ),
                      deltacp = 0):
        if (self.oscil_params['dm21'] != dm21) or (self.oscil_params['dm31'] != dm31) or \
                (self.oscil_params['theta23'] != theta23) or (self.oscil_params['theta13'] != theta13) or \
                (self.oscil_params['theta12'] != theta12) or (self.oscil_params['deltacp'] != deltacp):
            self.changed = True
            
            self.oscil_params['dm21'] = dm21
            self.oscil_params['dm31'] = dm31
            self.oscil_params['theta23'] = theta23
            self.oscil_params['theta12'] = theta12
            self.oscil_params['theta13'] = theta13
            self.oscil_params['deltacp'] = deltacp            
            


    #-----------------------------------------------------------------------
    # Get the indicies for each MC set ONCE and save them with that set 
    # somehow. That's how we'll look up the probabilities
    #-----------------------------------------------------------------------
    def get_indicies(self, energy, coszenith, pdg_encoding):
        ptype = np.abs(pdg_encoding[0])
        
        # Numpy is kind of stupid here.
        # Need to do one index at a time.
        index_e = np.digitize( energy, bins = self.ebins)-1
        index_cz = np.digitize( coszenith, bins = self.czbins)-1
        index_p = np.digitize( pdg_encoding, self.nubar*ptype)-1
        
        return index_e, index_cz, index_p

    #-----------------------------------------------------------------------
    # Get the actual probabilities by interfacing with prob3
    # Only recalculate them if needed
    # Note: if you're using indicies, specify a scalar (integer) pdg!
    #-----------------------------------------------------------------------
    def get_probabilities(self, e=[], cz=[], pdg=[], indicies=None, params=None):
        if params: self.setParameters(**params)

        if self.changed:
            for particle in [12, 14, 16]:
                pmap = getProb3(self.emap.flatten(), 
                                self.czmap.flatten(), 
                                self.nu_nubar_map.flatten() * particle,
                                self.oscil_params)
                
                self.pmap_from_e[particle] = np.reshape(pmap[:,0],
                                                        self.emap.shape)
                self.pmap_from_mu[particle] = np.reshape(pmap[:,1],
                                                         self.emap.shape)

                # Spline it for smoothing
                for nubar in [-1,1]:
                    if nubar == -1: i = 0
                    else: i = 1

                    self.pspline[nubar*12][nubar*particle] = RectBivariateSpline( self.ebins,
                                                                                  self.czbins,
                                                                                  self.pmap_from_e[particle][:,:,i])
                    self.pspline[nubar*14][nubar*particle] = RectBivariateSpline( self.ebins,
                                                                                  self.czbins,
                                                                                  self.pmap_from_mu[particle][:,:,i])

            self.changed = False

        # Okay! I have the probabilities for each bin.
        # Now evaluate and return the probabilities
        if indicies:
            return np.array([self.pmap_from_e[pdg][indicies], 
                             self.pmap_from_mu[pdg][indicies]]).T

        # Otherwise use the splines for fanciness
        else:
            nubar = pdg < 0
            
            ptype = np.abs(pdg)[0]
            probs_e, probs_mu = np.ones_like(e), np.ones_like(e)
  
            probs_e[nubar] = self.pspline[-12][-1*ptype].ev( e[nubar], cz[nubar], )
            probs_e[~nubar] = self.pspline[12][ptype].ev( e[~nubar], cz[~nubar], )

            probs_mu[nubar] = self.pspline[-14][-1*ptype].ev( e[nubar], cz[nubar], )
            probs_mu[~nubar] = self.pspline[14][ptype].ev( e[~nubar], cz[~nubar], )
            
            return np.array( [probs_e, probs_mu] ).T
        

##########################################################################
# Test the class and function
##########################################################################
if __name__ == "__main__":
    energies = np.logspace(0, 2, 500000)
    coszens = -1#np.linspace(-1, 0, 1)
    pdg = -12

    energies, coszens, pdg = np.meshgrid(energies, coszens, pdg, indexing='ij')
    energies = energies.flatten()
    coszens = coszens.flatten()
    pdg = pdg.flatten()
    
    
    params = {'dm21' : 7.49e-5,
              'dm31' : 2.526e-3,
              'theta23' : np.arcsin( np.sqrt(0.440) ),
              'theta12' : np.arcsin( np.sqrt(0.308) ),
              'theta13' : np.arcsin( np.sqrt(0.002163) ),
              'deltacp' : 0}


    import time

    probcalc = binned_prob3()

    # Get the normal raw stuff from prob3 directly
    start = time.time()
    p3probs = getProb3(energies, coszens, pdg, params)
    ptime = time.time() - start
    print '500000 prob3 evaluations took ', str(ptime), 'seconds'

    # test using the indicies
    start = time.time()
    indicies = probcalc.get_indicies(energies, coszens, pdg)
    myprobs = probcalc.get_probabilities(pdg=np.abs(pdg[0]), indicies=indicies, params=params)
    itime = time.time() - start
    print '500000 table lookups took ', str(itime), 'seconds'
    print '\t Mean deviation:', np.std( (p3probs-myprobs)[energies>5] )
    print '\t Max deviation:', np.max( (p3probs-myprobs)[energies>5] )

    # And splines
    probcalc.changed = True
    start = time.time()
    spprobs = probcalc.get_probabilities(e=energies, cz=coszens, pdg=pdg, params=params)
    stime = time.time() - start
    print '500000 spline lookups took ', str(stime), 'seconds'
    print '\t Mean deviation:', np.std( (p3probs-spprobs)[energies>5] )
    print '\t Max deviation:', np.max( (p3probs-spprobs)[energies>5] )
    

    # Make a plot to show the various methods in e->mu and mu->mu probabilities
    from matplotlib import pyplot
    pyplot.figure()
    pyplot.plot( energies, myprobs[:,0], label = 'from indicies, e->mu', linewidth=2, alpha=0.5, linestyle='dotted')
    pyplot.plot( energies, myprobs[:,1], label = 'from indicies, mu->mu', linewidth=2, alpha=0.5, linestyle='dotted')
    pyplot.plot( energies, spprobs[:,0], label = 'from splines, e->mu', linewidth=2, alpha=0.5, linestyle='dashed')
    pyplot.plot( energies, spprobs[:,1], label = 'from splines, mu->mu', linewidth=2, alpha=0.5, linestyle='dashed')
    pyplot.plot( energies, p3probs[:,0], label = 'from full p3, e->mu', linewidth=2, alpha=0.5)
    pyplot.plot( energies, p3probs[:,1], label = 'from full p3, mu->mu', linewidth=2, alpha=0.5)
    pyplot.xscale('log')
    pyplot.ylim(0., 1)
    pyplot.grid()
    pyplot.legend()
    pyplot.savefig('prob3check.pdf')

    print 'Complete!'
    
