#bk#!/usr/bin/env python                                                                                
from os import uname
import argparse,os
import numpy 
from icecube import icetray, dataio, dataclasses, linefit, phys_services
from icecube import gulliver, gulliver_modules, spline_reco, StartingTrackVeto
from icecube.DeepCore_Filter import DOMS
from I3Tray import I3Units
from I3Tray import *
#from weighting import weighter
import STV_utilities as Utils
import STV_modules as Mods
import STV_cuts as Cuts
print uname()

dlist = DOMS.DOMS("IC86")
parser = argparse.ArgumentParser()

icetray.I3Logger.global_logger = icetray.I3PrintfLogger()
icetray.logging.console()
parser.add_argument("-i","--infiles",
                    dest="infiles",
                    type=str,
                    default=[],
                    nargs="+",
                    help="[I]nfiles with frames")

parser.add_argument("-o","--outfile",
                    dest="outfile",
                    type=str,
                    default="",
                    help="base name for [o]utfiles")

parser.add_argument("-v","--verbose",
                    dest="verbose",
                    type=int,
                    default=0,
                    help="verbose? 1=True 0=False")

parser.add_argument("-p","--pulses",
                    dest="pulses",
                    type=str,
                    default="SplitInIcePulses",
                    help="Name of pulse series to use.")

parser.add_argument("-sp","--srt_pulses",
                    dest="srt_pulses",
                    type=str,
                    default="SRTInIcePulses",
                    help="Name of pulse series to use.")

parser.add_argument("-t","--pm_thrd",
                    dest="pm_thrd",type=float,
                    default=1,
                    help="veto threshold below which events should be kept, initial veto calc cut")

parser.add_argument("-f","--seed_fit",
                    dest="seed_fit",
                    type=str,
                    default="SPEFit2",
                    help="Fit to seed the detector scans with")

parser.add_argument("-k","--keep_gcd",
                    dest="keep_gcd",
                    type=int,
                    default=0, 
                    help="0: Don't save GCD with file 1: Do save GCD with file")

parser.add_argument('--sk', '--skip', 
                    dest='skip', 
                    type=int,
                    default=0, 
                    help='Number of events to skip')

parser.add_argument('-id','--evt_id', 
                    dest='event_id', 
                    type=int,
                    default=-1, 
                    help='Only run this on the event \
                    with the specified event id; \
                    if not specified, default is -1, \
                    run on all input events')

parser.add_argument('--ne', '--n_events', 
                    dest='n_events', 
                    type=int,
                    default=-1, 
                    help='Number of events to run.')

parser.add_argument('--nf', '--n_files', 
                    dest='n_files', 
                    type=int,
                    default=6, 
                    help='Number of files for weighting.')

args = parser.parse_args()

infiles=args.infiles
pulses=args.pulses
srt_pulses=args.srt_pulses
outfile=args.outfile
verbose=args.verbose
pm_thrd=args.pm_thrd
seed_fit=args.seed_fit
keep_gcd=args.keep_gcd
n_skip=args.skip
event_id=args.event_id
n_events=args.n_events
n_files=args.n_files

############################Parameters and Global vars ##############################
###BadOMs##
badOMs=[]

###ProcessEvent Module
#Counter for frames
Count = 0 

###FreeRecoParameters
#Reconstruction for first simple STV run 
d_step_f = 1
ang_step_f = 10
llh ="MPE"
safe_fit ="SafeFit"
log_name_f = "_Free_{0}_{1}".format(d_step_f,ang_step_f) 
free_fit_name = "Spline{0}_{1}{2}".format(llh,seed_fit,log_name_f)  
safe_fit_name = "Spline{0}_{1}{2}".format(llh,safe_fit,log_name_f)  
#print "Free fit name = ", free_fit_name
#print "Safe fit name = ", safe_fit_name

###Primary##
prim = "Primary"
d_step_p = 0
ang_step_p = 0
log_name_p ="_{0}_{1}".format(d_step_p,ang_step_p) 
prim_fit_name = "Spline{0}_{1}{2}".format(llh,prim,log_name_p)  

#First STV on FreeReco
n_segments_free = 1
dtype_free = "cherdat"
min_cad_dist = 150

###VetoFitRecoParameters
#Reconstruction for each veto hit
d_step_vf = 1
ang_step_vf = 1
log_name_vf ="_{0}_{1}".format(d_step_vf,ang_step_vf) 
N = 5

###GetTrackNames
#TrackNames = ["VetoFit_1211","VetoFit_6248","VetoFit_6850"]
TrackNames = Utils.GetTrackNames()
CorTrackNames = []
CorTrackNames = Utils.GetCorridorTrackNames(Tracks = CorTrackNames)
FitNames = Utils.GetFitNames(TrackNames, Llh =llh, DStep = d_step_vf, AngStep = ang_step_vf)
SafeFitsSearch = ["MPEFit","SPEFit2","LineFit","SafeLineFit"]
inf_muon_service = Utils.GetPhotonicsService(service_type="inf_muon")
#weight_calc = weighter(mctype = 'genie', nfiles = n_files) # or mctype='corsika'
#print TrackNames, FitNames
values = []
tray = I3Tray()
################################END####################################################

def ProcessEvent(frame,NSkip,NEvents): #Has to be here, beacause it uses global var count
    #    print frame["I3EventHeader"]  #Needed to be used with aci scripts
    global Count
    Count = Count + 1
    if Count <= NSkip and not (NSkip==0 and NEvents==-1):
#        print('Skipping event %d' % (Count) )
        return False
    if NEvents > 0:
        if Count <= NSkip+NEvents:
            print('Running on event %d' % (Count))
        return Count <= (NSkip+NEvents)
    print('Running on event %d' % (Count) )
    return True



def GetPmiss(frame, FitName):
    global values
    if frame.Has(FitName):
        for k in frame.keys():
            if ("prob_obs_0s" in k) and (FitName in k):
                pm = frame[k].value
                values.append(pm)

                if len(values)%1000 == 0:
                    print len(values)/1000

    return

def Weight(frame):
    weight_calc = weighter(mctype = 'genie', nfiles = 2) # or mctype='corsika'
    w = weight_calc(frame)
    frame["Weight"] = dataclasses.I3Double(w[0])
    print "weight = ", w[0]
    return

################################START################################################          

tray.Add("I3Reader","reader", FilenameList=infiles)
#For running only a certain numer of events and skipping
tray.Add(ProcessEvent, "ProcessEvent",  
         NSkip=n_skip, 
         NEvents=n_events)

tray.Add(Utils.GetBadDOMList,"BadDOMList",Streams=[icetray.I3Frame.DetectorStatus])
tray.AddModule("I3LineFit", "JustInCase",
              Name="SafeLineFit",
              InputRecoPulses="SplitInIcePulses",
              AmpWeightPower=1.)

tray.AddModule(Utils.SetSafeFit,"SetSafeFit", Fits = SafeFitsSearch, SafeFitName = safe_fit)
tray.AddModule("I3OMSelection<I3RecoPulseSeries>", 'selectSRTDCFidDOMs', #Select Fid DOMs
               OmittedKeys= dlist.DeepCoreFiducialDOMs,
               SelectInverse = True,
               InputResponse = srt_pulses,
               OutputResponse = 'SRTPulsesFid',
               OutputOMSelection = 'Selection_DCFidSRT',
               )

tray.AddModule("I3OMSelection<I3RecoPulseSeries>", 'selectSRTDCVetoDOMs', #Select Veto DOMs
               OmittedKeys= dlist.DeepCoreFiducialDOMs,
               SelectInverse = False,
               InputResponse = srt_pulses,
               OutputResponse = 'SRTPulsesVeto',
               OutputOMSelection = 'Selection_DCVetoSRT',
               )

tray.AddModule("I3OMSelection<I3RecoPulseSeries>", 'selectDCFidDOMs', #Select Fid DOMs
               OmittedKeys= dlist.DeepCoreFiducialDOMs,
               SelectInverse = True,
               InputResponse = pulses,
               OutputResponse = 'PulsesFid',
               OutputOMSelection = 'Selection_DCFid',
               )

tray.AddModule("I3OMSelection<I3RecoPulseSeries>", 'selectDCVetoDOMs', #Select Veto DOMs
               OmittedKeys= dlist.DeepCoreFiducialDOMs,
               SelectInverse = False,
               InputResponse = pulses,
               OutputResponse = 'PulsesVeto',
               OutputOMSelection = 'Selection_DCVeto',
               )


tray.AddModule(Cuts.CalculateVars, "CalcVars", Pulses = srt_pulses, PulsesFid = 'SRTPulsesFid', PulsesVeto = 'SRTPulsesVeto')
tray.AddModule(Cuts.PreCut, "PreCuts")
tray.AddModule(Mods.CoGMedIC, "CoG", PulsesFid = 'SRTPulsesFid')
#tray.AddModule(Weight,"Weight")
#tray.AddModule(Cuts.Primary,"Energy",FitName=safe_fit)

###########SimpleRecoSTV############
tray.AddSegment(Mods.DoVetoFits,"DoSplineRecoFit_Free", #DO a Fit on a seed_fit
                Pulses="SRTPulsesFid",
                FitNames=[seed_fit],
                Llh=llh,
                LogName=log_name_f,
                Spline=inf_muon_service,
                AngStep=ang_step_f,
                DistStep=d_step_f,
                If = lambda frame: frame.Has(seed_fit))


# ##Corridors
# tray.Add(Mods.MakeCorridorFits,"MakeCorridorFits", 
#          CoGFit = free_fit_name,
#          SafeFit = safe_fit)

# tray.AddSegment(Mods.DoTrackHits,"DoTH_Corrs", 
#                 Pulses=pulses,
#                 FitNames=CorTrackNames,
#                 NSeg=n_segments_free,
#                 Percent=0.01,
#                 Spline=inf_muon_service,
#                 MinCADist=min_cad_dist)

# tray.Add(Utils.CleanTH,"CleanTHCorFits", 
#          Pulses=pulses, 
#          FitNames=CorTrackNames) 

# tray.AddSegment(Mods.EvalLLH,"EvalLLHCorrs", 
#          Llh=llh, 
#          FitNames=CorTrackNames, 
#          Pulses = pulses,
#          Spline = inf_muon_service)

# tray.Add(Mods.SelectLLH,"SelectLLHCorrs", 
#          TrackNames=CorTrackNames, 
#          N = N)


# #############Search###############

tray.Add(Mods.MakeVetoFits,"MakeVetoFits", 
         CoGFit = free_fit_name,
         SafeFit = safe_fit,
         PulsesVeto = 'PulsesVeto')

tray.AddSegment(Mods.DoTrackHits,"DoTH_Veto", 
                Pulses=pulses,
                FitNames=TrackNames,
                NSeg=n_segments_free,
                Percent=0.01,
                Spline=inf_muon_service,
                MinCADist=min_cad_dist)



# tray.Add(Utils.CleanTH,"CleanTHVetoFits", 
#          Pulses=pulses, 
#          FitNames=TrackNames) 

tray.Add(Mods.MakeVetoPulses,"MakeVetoPulses",
         PulsesVeto = 'PulsesVeto',
         PulsesFid = 'PulsesFid')

tray.AddSegment(Mods.EvalLLHInit,"EvalLLHInit", 
         Llh=llh, 
         FitNames = TrackNames, 
         Spline = inf_muon_service)


tray.AddSegment(Mods.DoVetoPulseFits, "DoVetoPulseFits",
                FitNames = TrackNames,
                Llh = llh,
                LogName = log_name_vf,
                Spline = inf_muon_service,
                AngStep = ang_step_vf,
                DistStep = d_step_vf)

tray.AddSegment(Mods.EvalLLHFin,"EvalLLHFin", 
         Llh=llh, 
         FitNames = FitNames, 
         Spline = inf_muon_service)


# tray.Add(Utils.CleanReco,"CleanRecoVetoFits", 
#          Llh=llh, 
#          FitNames = FitNames, 
#          AngStep = ang_step_vf,
#          DistStep = d_step_vf)

tray.AddSegment(Mods.EvalLLH,"EvalLLH", 
         Llh=llh, 
         FitNames = FitNames, 
         Pulses = pulses,
         Spline = inf_muon_service)


# #tray.Add(Utils.PrintLLH,"PrintLLHVetofits2", FitNames = FitNames)
tray.Add(Mods.SelectLLH,"SelectLLH", 
         TrackNames=TrackNames, 
         N = N)
# #Find P_miss
tray.AddSegment(Mods.DoSTV,"DoSTV_VetoFits",
                Pulses=pulses,
                FitNames=FitNames,
                NSeg=n_segments_free,
                PmCut=pm_thrd,
                Spline=inf_muon_service,
                MinCADist=min_cad_dist,
                DistType=dtype_free)

tray.Add(Utils.CleanSTV,"CleanSTVVetoFits", 
         Pulses=pulses, 
         FitNames = FitNames+CorTrackNames)

tray.Add(Utils.PrintLLH,"PrintLLHVetofits", FitNames = TrackNames)

##########Writing I3 files#############
if len(outfile)==0:
    outfile="test"
if keep_gcd:
    tray.AddModule("I3Writer", "writer", Filename=outfile+".i3.bz2",
                   Streams=[icetray.I3Frame.Geometry,
                            icetray.I3Frame.Calibration,
                            icetray.I3Frame.DetectorStatus,
                            icetray.I3Frame.DAQ,
                            icetray.I3Frame.Physics],
                   DropOrphanStreams=[icetray.I3Frame.DAQ])
else:
    tray.AddModule("I3Writer", "writer", Filename=outfile+".i3.bz2",
                   Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
                   DropOrphanStreams=[icetray.I3Frame.Geometry,
                                      icetray.I3Frame.Calibration,
                                      icetray.I3Frame.DetectorStatus,
                                      icetray.I3Frame.DAQ])

tray.AddModule("TrashCan", "thecan")
#tray.Execute()
tray.Execute(4+2*(n_skip+n_events))
tray.Finish()

