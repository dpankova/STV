#!/usr/bin/python
from __future__ import division
import icecube
from I3Tray import *
from operator import itemgetter
import numpy 
import glob
import sys
from icecube import icetray, phys_services, dataclasses, dataio, photonics_service, TrackHits
import pickle
import copy 


file_list = []
name_f = "TrackHitsResult"
data_file = "/gpfs/group/dfc13/default/dasha/mlarson/L2/nutau/genie_ic.16640.000000.i3.gz"
gcd_file = "/gpfs/group/dfc13/default/dasha/mlarson/L2/GeoCalibDetectorStatus_2013.56429_V1_Modified.i3.gz"
file_list.append(gcd_file)
# for filename in glob.glob(data_file):
file_list.append(data_file)

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


inf_muon_service = GetPhotonicsService(service_type="inf_muon")

#Clean up the debug from the frame and book max compatible hits
def CleanTH(frame, Pulses, TrackNames):
    lists = []
    fitnames = []
    for fitname in TrackNames:
        if frame.Has(fitname):
            p = 0
            q = 0
            for k in frame.keys():
                if ("TrackHits_{0}_{1}".format(fitname,Pulses) in k) and ("coincObsQsList" in k):
                    Ps = frame[k]
                    for om, charges in Ps:
                        if not len(charges) == 0:
                            p = p + 1
                            q = q + sum(charges)
                    del frame[k]

                elif ("TrackHits_{0}_{1}".format(fitname,Pulses) in k):
                    del frame[k]

            lists.append([p,q,fitname])
            
    ps = sorted(lists, key=itemgetter(0), reverse=True)
    qs = sorted(lists, key=itemgetter(1), reverse=True)
    frame["TrackHits_MaxCompHits"]=dataclasses.I3Double(ps[0][0])
    frame["TrackHits_MaxCompCharge"]=dataclasses.I3Double(qs[0][1])
    print "TrackHits_MaxCompHits = {0:.3f}".format(ps[0][0])
    print "TrackHits_MaxCompCharge = {0:.3f}".format(qs[0][1])
    
def NSegmentVector(frame,FitName,N=1):
    #Make Segments out of tracks, Required for STV aan TH (made by K.Jero)
    if frame.Has(FitName):
        if N%2==0:
            print "n=",N,"is even! Change this!"
            sys.exit(910)
        try:
            basep=copy.deepcopy(frame[FitName])
        except:
            return True
        basep.shape = basep.shape.InfiniteTrack
        ##shift to closest approach to 0,0,0
        origin_cap = phys_services.I3Calculator.closest_approach_position(
            basep,dataclasses.I3Position(0,0,0))
        basep_shift_d=numpy.sign(origin_cap.z - basep.pos.z)*numpy.sign(
            basep.dir.z)*(origin_cap-basep.pos).magnitude
        basep_shift_pos=basep.pos+(basep.dir*basep_shift_d)
        basep_shift_t=basep_shift_d/basep.speed
        basep.pos=basep_shift_pos
        basep.time=basep.time+basep_shift_t
        segments=[]
        segment_length=1950./N
        for idx in range(N):
            dshift=segment_length*(idx-((N-1)/2.))
            particle=dataclasses.I3Particle()
            particle.time=basep.time+(dshift/basep.speed)
            particle.pos=basep.pos+(basep.dir*dshift)
            particle.dir=basep.dir
            particle.energy=0.01
            if N==1:
                particle.shape=particle.shape.InfiniteTrack
                particle.length=0
            else:
                particle.shape=particle.shape.ContainedTrack
                particle.length=segment_length
            segments.append(particle)

        del frame[FitName+"_"+str(N)+"_segments"]
        frame[FitName+"_"+str(N)+"_segments"]=dataclasses.I3VectorI3Particle(segments)
#        print "Put", FitName+"_"+str(N)+"_segments", "in the frame"
        del segments

def DoTrackHits(tray, name, Pulses, FitNames, NSeg, Percent, Spline, MinCADist):
    #Run TrackHits Module, looking for compatible Hits
    for fitname in FitNames:
        tray.AddModule(NSegmentVector,"NSegmentVectorTH_"+fitname+"_"+str(NSeg),
                       FitName=fitname,
                       N=NSeg)
        tray.Add("TrackHits","TH_"+fitname+"_"+str(NSeg),
                 Pulses=Pulses,
                 Photonics_Service=Spline,
                 Percent=Percent,
                 Fit=fitname,
                 Particle_Segments=fitname+"_"+str(NSeg)+"_segments",
                 Min_CAD_Dist=MinCADist)
    tray.Add(CleanTH,"CleanTH",
             TrackNames=FitNames,
             Pulses=Pulses)
        
#@icetray.traysegment
tray = I3Tray()
tray.AddModule("I3Reader","reader", FilenameList = file_list)
tray.AddSegment(DoTrackHits,"DoTH",
                Pulses="SplitInIcePulses",
                FitNames=["LineFit", "MPEFit"],
                NSeg=1,#Number of segments to split the track in, odd, 
                Percent=0.01,#Percent of expected Probability max, defines size of time window where we look for pulses 
                Spline=inf_muon_service, #load Tables
                MinCADist=150) #Look for hits inside this radius

tray.AddModule('I3Writer', 'writer', Filename= name_f+'.i3.bz2')
tray.AddModule('TrashCan','thecan')
tray.Execute(10)
tray.Finish()
