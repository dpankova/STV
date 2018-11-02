import os
import numpy
import pickle
import copy
from icecube import photonics_service, dataclasses

def GetBadDOMList(frame):
    global BadOMs
    BadOMs=frame["BadDomsList"]
    BadOMs.extend(frame["BadDomsListSLC"])

def CheckSeedFit(frame, seed_fit=""):
    if frame.Has(seed_fit):
#        fit = frame[seed_fit]
#        if fit.fit_status == dataclasses.I3Particle.OK:
        return True

    print "Seed is bad"
    return False

def SetSafeFit(frame, Fits, SafeFitName):
    # you can order fits
    for name in Fits:
        if frame.Has(name):
            fit = copy.deepcopy(frame[name])
            if fit.fit_status == dataclasses.I3Particle.OK: 
                frame[SafeFitName] = fit
                return True
    
    print "No Safe Fit"
    return    

def GetCorridorTrackNames(Tracks):
    inputs = open('Corridors.pkl', 'rb')
    data = pickle.load(inputs)
    for cor in data:
        for trk in data[cor][1]:
            Tracks.append(trk[0])
    return Tracks

def GetBestFineFitsNames():
    Tracks = []
    for i in range(0,5):
        name = "BestFineFits_"+str(i) 
        Tracks.append(name)
     
    return Tracks

def GetBestCoarseFitsNames():
    Tracks = []
    for i in range(0,5):
        name = "BestCoarseFits_"+str(i) 
        Tracks.append(name)
     
    return Tracks

    
# def GetTrackNames():
# #    tracknames=["LineFit_DC","SPEFit2","MPEFit"]
#     tracknames=[]
#     ic_om = range(1,61)
#     ic_st = range(1,79)
#     ic_om_names = ['VetoFit_{0:02d}{1:02d}'.format(st,om) for om in ic_om for st in ic_st]
#     dc_st = [26,27,35,36,37,45,46]
#     dc_om = range(38,61)
#     dc_om_names = ['VetoFit_{0:02d}{1:02d}'.format(st,om) for om in dc_om for st in dc_st]
#     set_dc = set(dc_om_names)
#     for i in ic_om_names:
#         if i not in set_dc:
#             tracknames.append(i)

#     return tracknames

def GetTrackNames():
#    tracknames=["LineFit_DC","SPEFit2","MPEFit"]
    tracknames=[]
    lim = 70
    for i in range(0,lim):
        name = 'VetoFit_{0:05d}'.format(i)
        tracknames.append(name)

    return tracknames

# def GetPulsesNames():
# #    tracknames=["LineFit_DC","SPEFit2","MPEFit"]
#     tracknames=[]
#     ic_om = range(1,61)
#     ic_st = range(1,79)
#     ic_om_names = ['VetoFit_{0:02d}{1:02d}_Pulses'.format(st,om) for om in ic_om for st in ic_st]
#     dc_st = [26,27,35,36,37,45,46]
#     dc_om = range(38,61)
#     dc_om_names = ['VetoFit_{0:02d}{1:02d}_Pulses'.format(st,om) for om in dc_om for st in dc_st]
#     set_dc = set(dc_om_names)
#     for i in ic_om_names:
#         if i not in set_dc:
#             tracknames.append(i)

#     return tracknames

def GetFineTrackNames(NPositionsFid, NPositionsVeto):
    tracknames = []
    for i in range(0,5):
        for index in range(0,NPositionsFid*NPositionsVeto):
            tracknames.append('FineFit_{0:02d}_{1:04d}'.format(i,index))

    return tracknames

def GetFitNames(TrackNames, Llh, DStep, AngStep):
    fitnames = []
    for track in TrackNames:
        name = "Spline{0}_{1}_{2!s}_{3!s}".format(Llh, track, DStep, AngStep)
        fitnames.append(name)

    return fitnames

def CheckVetoFits(frame, FitNames):
    for fitname in FitNames:
        if frame.Has(fitname):
            print "HERE!", fitname

def CleanTracks(frame, FitNames):
    for fitname in FitNames:
        if frame.Has(fitname):
            del frame[fitname]

def CleanSTV(frame, Pulses, FitNames, CleanAll=False ):
    for fitname in FitNames:
        if frame.Has(fitname): 
            for k in frame.keys():
                if CleanAll == False:
                    if ("{0}_{1}".format(Pulses,fitname) in k) and not ("prob_obs_0s" in k) :
                        del frame[k]
                else:
                    if ("{0}_{1}".format(Pulses,fitname) in k):
                        del frame[k]

def CleanTH(frame, Pulses, FitNames, CleanAll=False ):
    for fitname in FitNames:
        for k in frame.keys():
            if CleanAll == False:
                if ("TrackHits_{0}_{1}".format(fitname,Pulses) in k)\
                        and not ("coincObsPsList" in k) and not ("coincObsQsList" in k)\
                        and not ("coincObsProbsList" in k) :
                    del frame[k]

                elif ("TrackHits_{0}_{1}".format(fitname,Pulses) in k):   
                    lists = frame[k]
                   # print lists
                    nz_OMs = []
                        
                    for om, val in lists: #Find non zero hits
                        if val:
                             nz_OMs.append(om)

                    if not nz_OMs: #Didn't find any non-zero hits
                        del frame[k]
            else:
                if ("TrackHits_{0}_{1}".format(fitname,Pulses) in k):
                    del frame[k]

def CleanSegments(frame, FitNames, N):
    for fitname in FitNames:
        for k in frame.keys():
            if ("{0}_{1}_segments".format(fitname, N) in k):
                del frame[k]

def CleanMisc(frame):
    remove = ['FarthestDistFitDOM', 'LatestTimeFitDOM', "CoGTime", "CoGPos", 'BestCoarseFitPos','BestCoarseFit', 'FineFitPosFid','FineFitPosVeto']
    for k in frame.keys():
        for i in remove:
            if (k == i):
                del frame[k]


def CleanFitParams(frame, FitNames):
    for fitname in FitNames:
        for k in frame.keys():
            if ("{0}{4}".format(fitname,"FitParams") in k):
                del frame[k]
                

def CleanReco(frame, Llh, FitNames, AngStep, DistStep, CleanAll = False):
    for fitname in FitNames:
        for k in frame.keys():
#            print fitname
            if ("Spline{0}_{1}_{2!s}_{3!s}{4}".format(Llh,fitname,AngStep,DistStep,"FitParams") in k):
                fit_params = frame[k]    
 #               print fit_params
                if CleanAll == False:
                    if numpy.isnan(fit_params.rlogl):
                        del frame[k]
                        del frame[k[:-len("FitParams")]]
                else:
                        del frame[k]
                        del frame[k[:-len("FitParams")]]

def CleanEval(frame, Llh, FitNames, AngStep, DistStep, CleanAll = False):
    for fitname in FitNames:
        for k in frame.keys():
            #            print fitname
            if ("LLHCalc{0}_Spline{0}_{1}_{2!s}_{3!s}".format(Llh,fitname,AngStep,DistStep) in k):
                fit_params = frame[k]    
#                print fit_params
                if CleanAll == False:
                    if numpy.isnan(fit_params.rlogl):
                        del frame[k]
                else:
                    del frame[k]
                    
                        
def PrintLLH(frame, FitNames):
    for fitname in FitNames:
        if frame.Has(fitname):
            logl = 0  
            pm = 0
            ps = 0
            for k in frame.keys():
               # if ("prob_obs_0s" in k) and (fitname in k):
                #    pm = frame[k].value
                if ("LLHCalcMPE" in k) and (fitname in k):
                    logl = frame[k].logl
                if ("prob_obs_0s" in k) and (fitname in k):
                    pm = frame[k].value
                if ("coincObsQsList" in k) and (fitname in k):
                    Ps = frame[k]
                    for om, value in Ps:
                        if value:
                            ps = ps + 1

            if pm == 0:        
                print "{0} LLH = {1:.3e}".format(fitname, logl)
            else:
                print "{0} LLH = {1:.3e} Pm = {2:.3e} Pulses = {3:f}".format(fitname, logl, pm, ps)
#            print "Getting out", fitname, frame[fitname]

def GetPhotonicsService(service_type="inf_muon"):
    table_base=""
    if os.path.isfile(os.path.expandvars("$I3_DATA/photon-tables/splines/ems_mie_z20_a10.%s.fits") % "abs"):
        table_base = os.path.expandvars("$I3_DATA/photon-tables/splines/ems_mie_z20_a10.%s.fits")
    elif os.path.isfile("splines/ems_mie_z20_a10.%s.fits" % "abs"):
        table_base = os.path.expandvars("splines/ems_mie_z20_a10.%s.fits")
    elif os.path.isfile("/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/ems_mie_z20_a10.%s.fits" % "abs"):
        table_base = os.path.expandvars("/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/ems_mie_z20_a10.%s.fits")
    elif os.path.isfile("/home/icecube/i3/data/generalized_starting_events/splines/ems_mie_z20_a10.%s.fits" % "abs"):
        table_base = os.path.expandvars("/home/icecube/i3/data/generalized_starting_events/splines/ems_mie_z20_a10.%s.fits")
    else:
        print "You don't have splines anywhere I can find. This will eventually raise an error, for now it semi-silently dies"
    if service_type=="cscd":
        cascade_service = photonics_service.I3PhotoSplineService(table_base % "abs", table_base % "prob", 0,maxRadius    = 600.0)
        return cascade_service
    elif service_type=="seg_muon":
        seg_muon_service = photonics_service.I3PhotoSplineService(
                           amplitudetable = os.path.join( os.path.expandvars("$I3_DATA/photon-tables/splines/") ,"ZeroLengthMieMuons_250_z20_a10.abs.fits"),  ## Amplitude tables 
                           timingtable = os.path.join( os.path.expandvars("$I3_DATA/photon-tables/splines/") ,"ZeroLengthMieMuons_250_z20_a10.prob.fits"),    ## Timing tables
                           timingSigma  = 0.0,
                           maxRadius    = 600.0)
        return seg_muon_service
    elif service_type=="inf_muon":
        inf_muon_service = photonics_service.I3PhotoSplineService(
                           amplitudetable = os.path.join( os.path.expandvars("$I3_DATA/photon-tables/splines/") ,"InfBareMu_mie_abs_z20a10.fits"),  ## Amplitude tables 
                           timingtable = os.path.join( os.path.expandvars("$I3_DATA/photon-tables/splines/") ,"InfBareMu_mie_prob_z20a10.fits"),    ## Timing tables
                           timingSigma  = 0.0,
                           maxRadius    = 600.0) 
        return inf_muon_service
    else:
        print "You didn't give me a spline service type I recognize. This will eventually raise an error, for now it semi-silently dies"
