##D.Pankova 02/28/2019
##Utilities for STV for Loe energies

import os
import numpy
import pickle
import copy
from icecube import photonics_service, dataclasses

def GetBadDOMList(frame):
    global BadOMs
    BadOMs=frame["BadDomsList"]
    BadOMs.extend(frame["BadDomsListSLC"])

#Get a Fit that is for sure good
def SetSafeFit(frame, Fits, SafeFitName):
    for name in Fits:
        if frame.Has(name):
            fit = copy.deepcopy(frame[name])
            if fit.fit_status == dataclasses.I3Particle.OK: 
                frame[SafeFitName] = fit
                return True
    
    print "No Safe Fit"
    return
    
#Get Names for Corridor Tracks, NOT USED
def GetCorridorTrackNames(Tracks):
    inputs = open('Corridors.pkl', 'rb')
    data = pickle.load(inputs)
    for cor in data:
        for trk in data[cor][1]:
            Tracks.append(trk[0])
    return Tracks

#Make Names for VetoTracks
def GetTrackNames():
    tracknames=[]
    lim = 70
    for i in range(0,lim):
        name = 'VetoFit_{0:05d}'.format(i)
        tracknames.append(name)
    return tracknames

#Make Names for VetoFits 
def GetFitNames(TrackNames, Llh):
    fitnames = []
    for track in TrackNames:
        name = "Spline{0}_{1}".format(Llh, track)
        fitnames.append(name)
    return fitnames

#Clean all but Pms and TH max
def CleanALL(frame, TrackNames):
    cln_keys = ["DC_CoG_Pos_STV", "DC_CoG_Time_STV", "Vars_Charge_Ratio",
                "Vars_Hits_Fid", "Vars_Pulses_Veto", "Vars_Vertex_Z",
                "Vars_Func_1", "Vars_Func_2", "SafeFit", "PulsesFid",
                "PulsesVeto", "SRTPulsesVeto", "STRPulsesFid"]
    
    for track in TrackNames:
        for k in frame.keys():
            if (track in k):
                del frame[k]

    for key in cln_keys:
        del frame[key]
                
#Clean out STV debug from the Frame 
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

#Clean out TreackHits Debug from the frame
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
                    nz_OMs = []                        
                    for om, val in lists: #Find non zero hits
                        if val:
                            nz_OMs.append(om)
                    if not nz_OMs: #Didn't find any non-zero hits
                        del frame[k]
            else:
                if ("TrackHits_{0}_{1}".format(fitname,Pulses) in k):
                    del frame[k]

#Clean out Reco from the frame                
def CleanReco(frame, Llh, FitNames, AngStep, DistStep, CleanAll = False):
    for fitname in FitNames:
        for k in frame.keys():
            if ("Spline{0}_{1}_{2!s}_{3!s}{4}".format(Llh,fitname,AngStep,DistStep,"FitParams") in k):
                fit_params = frame[k]    
                if CleanAll == False:
                    if numpy.isnan(fit_params.rlogl):
                        del frame[k]
                        del frame[k[:-len("FitParams")]]
                else:
                    del frame[k]
                    del frame[k[:-len("FitParams")]]

#Clean LLH evaluation from Frame
def CleanEval(frame, Llh, FitNames, AngStep, DistStep, CleanAll = False):
    for fitname in FitNames:
        for k in frame.keys():
            if ("LLHCalc{0}_Spline{0}_{1}_{2!s}_{3!s}".format(Llh,fitname,AngStep,DistStep) in k):
                fit_params = frame[k]    
                if CleanAll == False:
                    if numpy.isnan(fit_params.rlogl):
                        del frame[k]
                else:
                    del frame[k]
                    
#Print result                        
def PrintLLH(frame, FitNames):
    for fitname in FitNames:
        if frame.Has(fitname):
            logl = 0  
            pm = 0
            ps = 0
            for k in frame.keys():
                if ("LLHCalc" in k) and (fitname in k):
                    logl = frame[k].logl
                if ("prob_obs_0s" in k) and (fitname in k):
                    pm = frame[k].value
                if ("TrackHits" in k) and ("coincObsQsList" in k) and (fitname in k):
                    Qs = frame[k]
                    qs = []
                    for om,value in Qs:
                        if value:
                            qs.append(sum(value)) 
            if pm == 0:        
                print "{0} LLH = {1:.3e}".format(fitname, logl)
            else:
                print "{0} LLH = {1:.3e} Pm = {2:.3e} Pulses = {3:f}".format(fitname, logl, pm, sum(qs))

#Load Tables
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
