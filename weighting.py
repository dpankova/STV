#!/usr/bin/python
import icecube, numpy
from icecube import icetray, dataio, dataclasses
from icecube.dataclasses import I3Particle, I3Position
from icecube.icetray import I3Units

################################################################################
# wrapper to simplify weighting stuff
################################################################################
class weighter:
    def __init__(self, mctype='genie', nfiles=1, runs=[]):
        if 'muongun' in mctype.lower():
            self.weight = mg_weighter(nfiles)
        elif 'corsika' in mctype.lower():
            try: self.weight = corsika_weighter(runs, nfiles)
            except: 
                print "Got an error using the nice weighting code. Falling back on CorsikaWeightDict "
                self.weight = corsika_backup_weighter(runs, nfiles)

        elif 'genie' in mctype.lower():
            self.weight = neutrino_weighter(nfiles, False)
        elif 'nugen' in mctype.lower():
            self.weight = neutrino_weighter(nfiles, True)        
        elif 'noise' in mctype.lower():
            self.weight = noise_weighter(nfiles)
        elif 'data' in mctype.lower():
            self.weight = data_weighter(runs)

    def __call__(self, frame):
        return self.weight(frame)

    
################################################################################
# Weighting code for MuonGun
################################################################################
# class mg_weighter:
#     def __init__(self, nfiles):
#         model_name = "GaisserH4a_atmod12_SIBYLL"
        
#         # Set up the MuonGun weighting stuff
#         model = MuonGun.load_model(model_name)
#         model.flux.max_multiplicity = 1
#         outSurface = MuonGun.Cylinder(1600*I3Units.m, 800*I3Units.m)
#         inSurface = MuonGun.Cylinder(500*I3Units.m, 150*I3Units.m, I3Position(46.29,-34.88,-300))
#         spectrum = MuonGun.OffsetPowerLaw(5.0, 7e2, 160, 500)
#         scaling = MuonGun.ConstantSurfaceScalingFunction(inSurface)
#         generator = MuonGun.EnergyDependentSurfaceInjector(surface = outSurface,
#                                                            flux = model.flux,
#                                                            energy = spectrum, 
#                                                            radius = model.radius,
#                                                            scaling = scaling)
#         self.weighter = MuonGun.WeightCalculator(model, generator)
#         self.nfiles = nfiles
#         self.nevents = 1e5

#     def __call__(self, frame):
#         mctree = frame['I3MCTree']
#         muon = dataclasses.get_most_energetic_muon(mctree)
        
#         w = self.weighter([muon.pos.x], [muon.pos.y], [muon.pos.z],
#                         [muon.dir.zenith], [muon.dir.azimuth],
#                         [1.], [muon.energy], [0]) / ( self.nfiles * self.nevents)
#         return w        

################################################################################
# Weighting code for CORSIKA
################################################################################
from icecube import weighting
from icecube.weighting.fluxes import GaisserH4a
from icecube.weighting.weighting import from_simprod
class corsika_weighter:
    def __init__(self, runs, nfiles):
        self.generator = None
        self.flux = GaisserH4a()
        
        for i, run in enumerate(runs):
            current = from_simprod(run)
            if type(self.generator) == type(None): self.generator = current * nfiles[i]
            else: self.generator += current * nfiles[i]
            
        return

    def __call__(self, frame):
        mctree = frame['I3MCTree']
        primary = dataclasses.get_most_energetic_primary(mctree)
        
        w = self.flux(primary.energy, primary.pdg_encoding)/self.generator(primary.energy, primary.pdg_encoding)
        if numpy.isnan(w) or not numpy.isfinite(w): w = 0
        return w




################################################################################
# Weighting code for CORSIKA
################################################################################
from icecube import weighting
from icecube.weighting.fluxes import GaisserH4a
from icecube.weighting.weighting import from_simprod
class corsika_backup_weighter:
    def __init__(self, runs, nfiles):
        self.nfiles = nfiles
        return

    def __call__(self, frame):
        cwd = frame['CorsikaWeightMap']
        w = cwd['Weight'] / (cwd['TimeScale'] * self.nfiles[0])
        return w


################################################################################
# Weighting code for GENIE and NuGen
################################################################################
from icecube import NuFlux, prob3
from binned_prob3 import binned_prob3
class neutrino_weighter:
    def __init__(self, nfiles, isNuGen=False):
        self.lowe_flux_service = NuFlux.makeFlux("IPhonda2014_spl_solmin")
        self.highe_flux_service = NuFlux.makeFlux("honda2006")
        self.highe_flux_service.knee_reweighting_model = 'gaisserH3a_elbert'

        self.oscil_params = {"dm21": 7.49e-5,
                             "dm31": 2.526e-3,
                             "theta23": numpy.arcsin( numpy.sqrt(0.440) ),
                             "theta12": numpy.arcsin( numpy.sqrt(0.308) ),
                             "theta13": numpy.arcsin( numpy.sqrt(0.02163) ),
                             "deltacp": 289.*numpy.pi/180.,
                             }
        self.calc = binned_prob3()
        self.calc.setParameters( **self.oscil_params )

        self.nfiles = nfiles
        if isNuGen: self.gen_ratio=0.5
        else: self.gen_ratio=0.7

    def __call__(self, frame):
        mc_weights = frame['I3MCWeightDict']
        true_neutrino = dataclasses.get_most_energetic_neutrino(frame['I3MCTree'])
        
        true_energy = mc_weights['PrimaryNeutrinoEnergy']
        true_zenith = true_neutrino.dir.zenith
        true_azimuth = true_neutrino.dir.azimuth

        if true_neutrino.energy < 10000: flux_service = self.lowe_flux_service
        else: flux_service = self.highe_flux_service

        nue, numu = I3Particle.ParticleType.NuE, I3Particle.ParticleType.NuMu
        nu_nubar_genratio = self.gen_ratio
        if true_neutrino.pdg_encoding < 0:
            nue, numu = I3Particle.ParticleType.NuEBar, I3Particle.ParticleType.NuMuBar
            nu_nubar_genratio = 1-nu_nubar_genratio 

        flux_nue = flux_service.getFlux(nue, true_energy, numpy.cos(true_zenith))
        flux_numu = flux_service.getFlux(numu, true_energy, numpy.cos(true_zenith))

        p_osc = self.calc.get_probabilities(e=[true_energy,], cz=[numpy.cos(true_zenith),], 
                                            pdg=[true_neutrino.pdg_encoding,])

        one_weight = mc_weights['OneWeight']
        n_events = mc_weights['NEvents']
        norm = (1.0 / (n_events * self.nfiles * nu_nubar_genratio))

        w = norm * one_weight * (flux_nue * p_osc[:,0] + flux_numu * p_osc[:,1])
        return w

################################################################################
# weighting code for the accidental triggers
################################################################################
class noise_weighter:
    def __init__(self, nfiles=84300):
        self.nfiles = nfiles
        return

    def __call__(self, frame):
        return (1-2800*30e-6)/(self.nfiles * 1000.0 * (100e-3-20e-6))

################################################################################
# weighting code for data events
################################################################################
#from GoodRunSampler import getTime
class data_weighter:
    def __init__(self, runs):
		return
#        self.livetime = 0
#        for year in ['2012', '2013', '2014']:
#            self.livetime += getTime(runs, year, False)
#        return

#    def __call__(self, frame):
#        return 1.0/self.livetime
