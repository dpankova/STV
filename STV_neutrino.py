#!/usr/bin/env python                                                                                
from os import uname
import argparse,os
import numpy 
from icecube import icetray, dataio, dataclasses, linefit, phys_services
from icecube import gulliver, gulliver_modules, spline_reco, StartingTrackVeto
from icecube.DeepCore_Filter import DOMS
from I3Tray import I3Units
from I3Tray import *
from weighting import weighter
import STV_utilities as Utils
import STV_modules as Mods
import STV_cuts as Cuts
print uname()

dlist = DOMS.DOMS("IC86EDC")
parser = argparse.ArgumentParser()

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
                    help='Only run this on the event with the specified event id; if not specified, default is -1, run on all input events')

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
print "Free fit name = ", free_fit_name
print "Safe fit name = ", safe_fit_name

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
d_step_vf = 0
ang_step_vf = 0
log_name_vf ="_{0}_{1}".format(d_step_vf,ang_step_vf) 

###MPEFits for vertex study 
d_step_v = 1
ang_step_v = 1
vertex_seed = "MPEFit"
log_name_v = "_Vertex_{0}_{1}".format(d_step_v,ang_step_v) 
mpe_v_name = "Spline{0}_{1}{2}".format(llh,vertex_seed,log_name_v)  
log_name_v_srt = "_VertexSRT_{0}_{1}".format(d_step_v,ang_step_v) 
mpe_v_srt_name = "Spline{0}_{1}{2}".format(llh,vertex_seed, log_name_v_srt)  
log_name_v_fid = "_VertexFid_{0}_{1}".format(d_step_v,ang_step_v) 
mpe_v_fid_name = "Spline{0}_{1}{2}".format(llh,vertex_seed,log_name_v_fid)  
log_name_v_srt_fid = "_VertexSRTFid_{0}_{1}".format(d_step_v,ang_step_v) 
mpe_v_srt_fid_name = "Spline{0}_{1}{2}".format(llh,vertex_seed,log_name_v_srt_fid)  

###GetTrackNames
TrackNames = Utils.GetCoarseTrackNames()
TrackNames = Utils.GetCorridorTrackNames(Tracks = TrackNames)
FineTrackNames = Utils.GetFineTrackNames(NPositionsFid = 21, NPositionsVeto = 21)
FitNames = Utils.GetFitNames(TrackNames, Llh =llh, DStep = d_step_vf, AngStep = ang_step_vf)
FineFitNames = Utils.GetFitNames(FineTrackNames, Llh =llh, DStep = d_step_vf, AngStep = ang_step_vf)
BestFineFitsNames = Utils.GetBestFineFitsNames()
BestCoarseFitsNames = Utils.GetBestCoarseFitsNames()

SafeFitsSearch = ["MPEFit","SPEFit2","LineFit","SafeLineFit"]

inf_muon_service = Utils.GetPhotonicsService(service_type="inf_muon")

values = []
#weight_calc = weighter(mctype = 'genie', nfiles = n_files) # or mctype='corsika'

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

tray.Add(ProcessEvent, "ProcessEvent",  #For running only a certain numer of events and skipping
         NSkip=n_skip, 
         NEvents=n_events)

tray.Add(Utils.GetBadDOMList,"BadDOMList",Streams=[icetray.I3Frame.DetectorStatus])
tray.AddModule("I3LineFit", "JustInCase",
               Name="SafeLineFit",
               InputRecoPulses="SplitInIcePulses",
               AmpWeightPower=1.)

tray.AddModule(Utils.SetSafeFit,"SetSafeFit", Fits = SafeFitsSearch, SafeFitName = safe_fit)
#tray.Add(Utils.CheckSeedFit,"ChkSeedFit",seed_fit=seed_fit) #Check if seedfit exists
#tray.Add(Utils.CheckSeedFit,"ChkVertexFit",seed_fit=vertex_seed) #Check if seedfit exists
#tray.Add(Utils.PrintPm,"PrintPmSeed", FitNames = [seed_fit])
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


tray.AddModule(Cuts.CountAndCRatio, "PulseNCalc", PulsesFid = 'SRTPulsesFid', PulsesVeto = 'SRTPulsesVeto')
tray.AddModule(Cuts.PulsesCut, "PulseNCut", NumPFid = 5, NumPVeto = 100)
tray.AddModule(Mods.CoGMedIC, "CoG", PulsesFid = 'SRTPulsesFid')
tray.AddModule(Weight,"Weight")
tray.AddModule(Cuts.Primary,"Energy",FitName=safe_fit)
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

#Check if it worked
#tray.Add(Mods.RLogLCut,"RLoGLCut", fits_to_try=[free_fit_name])

# tray.AddSegment(Mods.DoSTV,"DoSTV_Free", #Run Starting TrackVeto on Seed_Fit
#                 Pulses=pulses,
#                 FitNames=[free_fit_name],
#                 NSeg=n_segments_free,
#                 PmCut=pm_thrd,
#                 Spline=inf_muon_service,
#                 MinCADist=min_cad_dist,
#                 DistType=dtype_free)

# tray.Add(Utils.CleanSTV,"CleanSTVFree",
#           Pulses=pulses,
#           FitNames=[seed_fit])

tray.Add(Utils.PrintPm,"PrintPmFree", FitNames = [free_fit_name])

#############SimpleRecoEND###############
################Vertex#################

# tray.AddSegment(Mods.DoVetoFits,"DoSplineRecoFit_Vertex1", #DO a Fit on a seed_fit
#                 Pulses=pulses,
#                 FitNames=[vertex_seed],
#                 Llh=llh,
#                 LogName=log_name_v,
#                 Spline=inf_muon_service,
#                 AngStep=ang_step_v,
#                 DistStep=d_step_v)

# tray.AddSegment(Mods.DoVetoFits,"DoSplineRecoFit_Vertex2", #DO a Fit on a seed_fit
#                 Pulses="SRTInIcePulses",
#                 FitNames=[vertex_seed],
#                 Llh=llh,
#                 LogName=log_name_v_srt,
#                 Spline=inf_muon_service,
#                 AngStep=ang_step_v,
#                 DistStep=d_step_v)

# tray.AddSegment(Mods.DoVetoFits,"DoSplineRecoFit_Vertex3", #DO a Fit on a seed_fit
#                 Pulses="PulsesFid",
#                 FitNames=[vertex_seed],
#                 Llh=llh,
#                 LogName=log_name_v_fid,
#                 Spline=inf_muon_service,
#                 AngStep=ang_step_v,
#                 DistStep=d_step_v)

# tray.AddSegment(Mods.DoVetoFits,"DoSplineRecoFit_Vertex4", #DO a Fit on a seed_fit
#                 Pulses="SRTPulsesFid",
#                 FitNames=[vertex_seed],
#                 Llh=llh,
#                 LogName=log_name_v_srt_fid,
#                 Spline=inf_muon_service,
#                 AngStep=ang_step_v,
#                 DistStep=d_step_v)
 
# tray.Add(Utils.PrintPm,"PrintPmVertex", FitNames = [mpe_v_name, mpe_v_srt_name, mpe_v_fid_name, mpe_v_srt_fid_name])

# ################Primary############################

tray.AddSegment(Mods.DoVetoFits,"DoSplineRecoFit_Prim", #DO a Fit on a seed_fit
                Pulses="SRTPulsesFid",
                FitNames=[prim],
                Llh=llh,
                LogName=log_name_p,
                Spline=inf_muon_service,
                AngStep=ang_step_p,
                DistStep=d_step_p)


tray.AddSegment(Mods.DoSTV,"DoSTV_Prim", 
                Pulses=pulses,
                FitNames=[prim+"Orig"],
                NSeg=n_segments_free,
                PmCut=pm_thrd,
                Spline=inf_muon_service,
                MinCADist=min_cad_dist,
                DistType=dtype_free)

tray.Add(Utils.CleanSTV,"CleanSTVPrimOrig",
         Pulses=pulses,
         FitNames=[prim+"Orig"])

tray.Add(Utils.PrintPm,"PrintPmPrim", FitNames = [prim+"Orig"])

tray.AddSegment(Mods.DoSTV,"DoSTV_PrimFit", 
                Pulses=pulses,
                FitNames=[prim_fit_name],
                NSeg=n_segments_free,
                PmCut=pm_thrd,
                Spline=inf_muon_service,
                MinCADist=min_cad_dist,
                DistType=dtype_free)

tray.Add(Utils.CleanSTV,"CleanSTVPrim",
         Pulses=pulses,
         FitNames=[prim])


tray.Add(Utils.PrintPm,"PrintPmPrimFit", FitNames = [prim_fit_name])

# #############CoarseGridSearch###############
tray.AddModule(Mods.CascadeDir, "CasDir", CoGFit = free_fit_name, SRTPulsesFid = 'SRTPulsesFid')

tray.Add(Mods.MakeVetoFits,"MakeVetoFits", 
         CoGFit = free_fit_name,
         SafeFit = safe_fit,
         PulsesVeto = 'PulsesVeto')

tray.Add(Mods.MakeCorridorFits,"MakeCorridorFits", 
         CoGFit = free_fit_name,
         SafeFit = safe_fit)

tray.AddSegment(Mods.DoVetoFits, "DoVetoFits",
                Pulses="SRTPulsesFid",
                FitNames = TrackNames,
                Llh = llh,
                LogName = log_name_vf,
                Spline = inf_muon_service,
                AngStep = ang_step_vf,
                DistStep = d_step_vf)

tray.Add(Utils.CleanReco,"CleanRecoVetoFits", 
         Llh=llh, 
         FitNames=TrackNames, 
         AngStep = ang_step_vf,
         DistStep = d_step_vf)

######CoarseFitSTV#########
# tray.AddSegment(Mods.DoSTV,"DoSTV_VetoFits",
#                 Pulses=pulses,
#                 FitNames=FitNames,
#                 NSeg=n_segments_free,
#                 PmCut=pm_thrd,
#                 Spline=inf_muon_service,
#                 MinCADist=min_cad_dist,
#                 DistType=dtype_free)

# tray.Add(Utils.CleanSTV,"CleanSTVVetoFits", 
#          Pulses=pulses, 
#          FitNames=FitNames)


#tray.Add(Utils.PrintPm,"PrintPmVetofits", FitNames = FitNames)
# tray.Add(Utils.PrintPm,"PrintPmVetof2its", FitNames = TrackNames)


# ########################FineGridSearch##########
tray.Add(Mods.GetBestCoarseFit,"FinefitPos", FitNames = FitNames)
tray.Add(Mods.MakeFineFitPos,"Coords", CoarseFits = BestCoarseFitsNames)
tray.Add(Mods.MakeFineFits,"Finefits", CoarseFits = BestCoarseFitsNames)

tray.AddSegment(Mods.DoVetoFits, "DoFineFits",
                Pulses="SRTPulsesFid",
                FitNames = FineTrackNames,
                Llh = llh,
                LogName = log_name_vf,
                Spline = inf_muon_service,
                AngStep = ang_step_vf,
                DistStep = d_step_vf)

tray.Add(Utils.CleanReco,"CleanRecoFineFits", 
         Llh=llh, 
         FitNames=FineTrackNames, 
         AngStep = ang_step_vf,
         DistStep = d_step_vf)

# #tray.Add(Utils.PrintPm,"PrintPmFinefits", FitNames = FineFitNames)

tray.Add(Mods.GetBestFineFit,"GetBestFineFit", FitNames=FineFitNames)

tray.AddSegment(Mods.DoSTV,"DoSTV_BestFineFit",
                Pulses=pulses,
                FitNames=BestFineFitsNames,
                NSeg=n_segments_free,
                PmCut=pm_thrd,
                Spline=inf_muon_service,
                MinCADist=min_cad_dist,
                DistType=dtype_free)

tray.Add(Utils.CleanSTV,"CleanSTVFineFit", 
         Pulses=pulses, 
         FitNames=BestFineFitsNames)

#tray.Add(Utils.PrintPm,"PrintPmFinefits", FitNames = FineFitNames)
tray.Add(Utils.PrintPm,"PrintPmBestCoarsefit", FitNames = BestCoarseFitsNames)
tray.Add(Utils.PrintPm,"PrintPmBestFinefit", FitNames = BestFineFitsNames)
#tray.Add(GetPmiss,"Values",  FitName = "BestFineFit")

# #########################CleanUp###############
tray.Add(Utils.CleanTracks,"CleanTracks", FitNames =TrackNames)
tray.Add(Utils.CleanTracks,"CleanFineTracks", FitNames =FineTrackNames)
tray.Add(Utils.CleanSegments,"CleanSegments", FitNames =FitNames, N = 1)
tray.Add(Utils.CleanSegments,"CleanFineSegments", FitNames =FineFitNames, N = 1)
tray.Add(Utils.CleanReco,"CleanRecoFitsAll", 
         Llh=llh, 
         FitNames=TrackNames, 
         AngStep = ang_step_vf,
         DistStep = d_step_vf,
         CleanAll = True)

tray.Add(Utils.CleanReco,"CleanRecoFineFitsAll", 
         Llh=llh, 
         FitNames=FineTrackNames, 
         AngStep = ang_step_vf,
         DistStep = d_step_vf,
         CleanAll = True)

# tray.Add(Utils.CleanMisc,"CleanMisc")

# ###########End#########


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

#numpy.save(outfile+".npy",numpy.array(values))                                               
