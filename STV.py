#!/usr/bin/env python                                                                                
#from os import uname
import argparse
import numpy 
from icecube import icetray, dataio, dataclasses, linefit, phys_services
from icecube import gulliver, gulliver_modules, spline_reco, StartingTrackVetoLE
from icecube.DeepCore_Filter import DOMS
from I3Tray import I3Units
from I3Tray import *
#from weighting import weighter

import STV_utilities as Utils
import STV_modules as Mods
import STV_cuts as Cuts


dlist = DOMS.DOMS("IC86")
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
outfile=args.outfile

n_skip=args.skip
event_id=args.event_id
n_events=args.n_events
n_files=args.n_files

#This modulse is needed here to be used with aci scripts
####################################
Count = 0 
def ProcessEvent(frame,NSkip,NEvents): 
    global Count
    Count = Count + 1
    if Count <= NSkip and not (NSkip==0 and NEvents==-1):
        #print('Skipping event %d' % (Count) )
        return False
    if NEvents > 0:
        if Count <= NSkip+NEvents:
            print('Running on event %d' % (Count))
        return Count <= (NSkip+NEvents)
    print('Running on event %d' % (Count) )
    return True



###################Parameters####################
#Those that are set manually:

#Pulses
PULSES="SplitInIcePulses"
SRTPULSES="SRTInIcePulses"

#LoadTables
TABLES = Utils.GetPhotonicsService(service_type="inf_muon")
#TABLES = Utils.GetPhotonicsService(service_type="seg_muon")

#For the fitter
LLH ="MPE" #Likelihood type: "SPE", "MPE", etc 
DSTEP = 1 #Distanse steps for minimizer
ASTEP = 1 #Angular steps ...

#There must be at least one fit with OK status to make tracks
#Doesn't matter how presise 
SAFEFITS = ["MPEFit","SPEFit2","LineFit"] #multiple candidates, incase some failed
SAFEFIT ="SafeFit" #Future name of the OK fit

#Vertex in DC, neede for Algorithm
#Future Name of the Vertex fit, which gives the vertex
VERTEXFIT = "Spline{0}_{1}".format(LLH,SAFEFIT)
#If Vertex fit fails, vertex = CoG in DC located here:
COGPOS = "DC_Cog_Pos_STV"
COGTIME = "DC_Cog_Time_STV"

#For STV
NSEG = 1 #Number of segments to split track into
DTYPE = "cherdat" #a way of calculating distance: cherdat, contribdat, caddat
RADIUS = 150 #Only look at DOMs closer than R to the track 

#TrackHitsCut threshold
#There is a plot in my talk
THCUT = 6

#Number of best tracks selected for Pm calculation
#The more the slower, but not really better
NTRACKS = 5

#List of track names for the algorithm
#Needs to be filled ahead of time, because modules 
#only interact with the frame, not each other
TRACKNAMES = Utils.GetTrackNames()
#Corresponding names of the fitted tracks 
#also needs to be created ahead of time
FITNAMES = Utils.GetFitNames(TRACKNAMES, Llh =LLH)

#For BadOMs modules
badOMs=[]

##############Algoritm##############              
def STV_LE_Algorithm(tray, name, TrackNames, FitNames,  
                     Pulses, PulsesFid, PulsesVeto,
                     VertexFit, CoGPos, CoGTime, SafeFit, 
                     THCutThd, LLH, AngStep, DistStep,
                     N, Tables, NSeg, R, DistType):
    
    #Make Tracks from Vertex to each Veto hit
    tray.Add(Mods.MakeVetoTracks,"MakeVetoTracks", 
             CoGFit = VertexFit,
             CoGPoSName = CoGPos,
             CoGTimeName = CoGTime,
             SafeFit = SafeFit,
             PulsesVeto = PulsesVeto)

    #Find compatible Pulses for each track
    tray.AddSegment(Mods.DoTrackHits,"RunTrackHits", 
                    Pulses=Pulses,
                    FitNames=TrackNames,
                    NSeg=NSeg,
                    Spline=Tables,
                    MinCADist=R)

    #Do a cut on max number of compatible hits
    #Comment out if you do it separately
    tray.Add(Cuts.THCut,"TrackHitsCut", 
             Pulses=Pulses, 
             FitNames=TrackNames,
             Threshold = THCutThd)
 
    #Make Pulses for Fitter = Compatible hits + DC Hits
    tray.Add(Mods.MakeVetoPulses,"MakeVetoPulses",
             PulsesVeto = PulsesVeto,
             PulsesFid = PulsesFid)

    #Run Fiiter on Tracks and corresponding Pulses
    tray.AddSegment(Mods.DoVetoPulseFits, "DoVetoPulseFits",
                    FitNames = TrackNames,
                    Llh = LLH,
                    Spline = Tables,
                    AngStep = AngStep,
                    DistStep = DistStep)

    #Eveluate LLH for Fits for comparision
    tray.AddSegment(Mods.EvalLLH,"EvalLLH", 
                    Llh=LLH, 
                    FitNames = FitNames, 
                    Pulses = Pulses,
                    Spline = Tables)

    #Select best tracks: N tracks with most comp Hits and lowest LLH
    tray.Add(Mods.SelectLLH,"SelectLLH", 
             TrackNames=TrackNames, 
             N = N)

    #Find Pm for best tracks
    tray.AddSegment(Mods.DoSTV,"DoSTV_VetoFits",
                Pulses=Pulses,
                FitNames=FitNames,
                NSeg=NSeg,
                Spline=Tables,
                MinCADist=R,
                DistType=DistType)

    #Print the result
    tray.Add(Utils.PrintLLH,"PrintLLHVetofits", 
             FitNames = TrackNames)
    
    #Book different types of Pm
    tray.Add(Mods.Pmiss,"PmissCalc", 
             Pulses=Pulses, 
             TrackNames=TrackNames)

    #Clean the frame (loots of extras stuff was added)
    tray.Add(Utils.CleanALL,"CleanALL", 
         TrackNames=TrackNames) 

#############EndAlgorithm##############################

#######Start###############
tray = I3Tray()
tray.Add("I3Reader","reader", FilenameList=infiles)

#For running only a certain numer of events and skipping
tray.Add(ProcessEvent, "ProcessEvent",  
         NSkip=n_skip, 
         NEvents=n_events)

#Get BadOMs
tray.Add(Utils.GetBadDOMList,"BadDOMList",Streams=[icetray.I3Frame.DetectorStatus])

####Calculate more inputs for the algorithm
#Split pulses into fiducial and veto pulses 
tray.AddModule("I3OMSelection<I3RecoPulseSeries>", 'selectSRTDCFidDOMs', #Select Fid DOMs
               OmittedKeys= dlist.DeepCoreFiducialDOMs,
               SelectInverse = True,
               InputResponse = SRTPULSES,
               OutputResponse = 'SRTPulsesFid',
               OutputOMSelection = 'Selection_DCFidSRT',
               )

tray.AddModule("I3OMSelection<I3RecoPulseSeries>", 'selectSRTDCVetoDOMs', #Select Veto DOMs
               OmittedKeys= dlist.DeepCoreFiducialDOMs,
               SelectInverse = False,
               InputResponse = SRTPULSES,
               OutputResponse = 'SRTPulsesVeto',
               OutputOMSelection = 'Selection_DCVetoSRT',
               )

tray.AddModule("I3OMSelection<I3RecoPulseSeries>", 'selectDCFidDOMs', #Select Fid DOMs
               OmittedKeys= dlist.DeepCoreFiducialDOMs,
               SelectInverse = True,
               InputResponse = PULSES,
               OutputResponse = 'PulsesFid',
               OutputOMSelection = 'Selection_DCFid',
               )

tray.AddModule("I3OMSelection<I3RecoPulseSeries>", 'selectDCVetoDOMs', #Select Veto DOMs
               OmittedKeys= dlist.DeepCoreFiducialDOMs,
               SelectInverse = False,
               InputResponse = PULSES,
               OutputResponse = 'PulsesVeto',
               OutputOMSelection = 'Selection_DCVeto',
               )

#Get a safe fit, what has OK status
tray.AddModule(Utils.SetSafeFit,"SetSafeFit", 
               Fits = SAFEFITS, 
               SafeFitName = SAFEFIT)

#Get a VertexFit, to get the vertex from
tray.AddSegment(Mods.DoVetoFits,"DoSplineRecoFit_Free", #DO a Fit on a seed_fit
                Pulses="SRTPulsesFid",
                FitNames=SAFEFIT,
                Llh=LLH,
                Spline=TABLES,
                AngStep=ASTEP,
                DistStep=DSTEP,
                If = lambda frame: frame.Has(SAFEFIT))

#If Vertexfit fails, get Vertex from CoG of DC pulses
tray.AddModule(Mods.CoGMedIC, "CoG", 
               PulsesFid = 'SRTPulsesFid',
               CoGPoSName = COGPOS,
               CoGTimeName = COGTIME)


#Calculate and use PreCuts
#Comment out if not using L2
tray.AddModule(Cuts.CalculateVars, "CalcVars", 
               SRTPulses = SRTPULSES, 
               SRTPulsesFid = 'SRTPulsesFid', 
               SRTPulsesVeto = 'SRTPulsesVeto',
               PulsesVeto = 'PulsesVeto')
tray.AddModule(Cuts.PreCut, "PreCuts")

#Run the Algoritm
tray.AddSegment(STV_LE_Algorithm, "TheAlgorithm", 
                TrackNames = TRACKNAMES, 
                FitNames = FITNAMES,  
                Pulses = PULSES, 
                PulsesFid = "PulsesFid", 
                PulsesVeto = "PulsesVeto",
                VertexFit = VERTEXFIT, 
                CoGPos = COGPOS, 
                CoGTime = COGTIME, 
                SafeFit  = SAFEFIT, 
                THCutThd = THCUT, 
                LLH = LLH, 
                AngStep = ASTEP, 
                DistStep = DSTEP,
                N = NTRACKS, 
                Tables= TABLES, 
                NSeg =NSEG, 
                R = RADIUS, 
                DistType = DTYPE)

if len(outfile)==0:
    outfile="test"

tray.AddModule("I3Writer", "writer", Filename=outfile+".i3.bz2",
               Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics],
               DropOrphanStreams=[icetray.I3Frame.Geometry,
                                  icetray.I3Frame.Calibration,
                                  icetray.I3Frame.DetectorStatus,
                                  icetray.I3Frame.DAQ])

tray.AddModule("TrashCan", "thecan")
tray.Execute(4+2*(n_skip+n_events))
tray.Finish()

