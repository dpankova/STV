#!/usr/bin/python
from __future__ import division
import icecube
from I3Tray import *
import numpy as np
import matplotlib
import matplotlib.path as mplPath
import glob
import sys
from math import sqrt, fabs, pi, tan, cos, sin
from icecube import icetray, dataclasses, dataio, phys_services, STTools
from icecube import photonics_service
from icecube import gulliver, gulliver_modules
from datetime import datetime
from icecube.icetray import I3Units
from argparse import ArgumentParser      
from icecube.DeepCore_Filter import DOMS
import pickle

#matplotlib.rcParams.update({'font.family': monospace})
matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams.update({'lines.linewidth' : 3})
dlist = DOMS.DOMS("IC86EDC")

#p1 = sys.argv[1]
#p2 = sys.argv[2]
   
file_list_n = []
file_list_c = []
file_list_d = []
gfile_n = '/gpfs/group/dfc13/default/dasha/mlarson/L2/GeoCalibDetectorStatus_2013.56429_V1_Modified.i3.gz'
gfile_c = '/gpfs/group/dfc13/default/dasha/mlarson/L2/corsika/GeoCalibDetectorStatus_2012.56063_V1.i3.gz'
gfile_d = '/gpfs/group/dfc13/default/dasha/mlarson/L2/data/Level2_IC86.2014_data_Run00125381_1004_0_116_GCD.i3.gz'
geofile = dataio.I3File(gfile_n)
file_list_n.append(gfile_n)
file_list_c.append(gfile_c)
file_list_d.append(gfile_d)

#data_file_1 = "/gpfs/group/dfc13/default/dasha/mlarson/pm_trial/nue/Pm_genie_ic.12640.00000{0}.part{1}*.i3.bz2".format(p1,p2)
#data_file_2 = "/gpfs/group/dfc13/default/dasha/mlarson/pm_trial/numu/Pm_genie_ic.14640.00000{0}.part{1}*.i3.bz2".format(p1,p2)
#data_file_3 = "/gpfs/group/dfc13/default/dasha/mlarson/pm_trial/nutau/Pm_genie_ic.16640.00000{0}.part{1}*.i3.bz2".format(p1,p2)
# data_file_1 = "/gpfs/group/dfc13/default/dasha/mlarson/pm_trial/nue/*.i3.bz2"
# data_file_2 = "/gpfs/group/dfc13/default/dasha/mlarson/pm_trial/numu/*.i3.bz2"
# data_file_3 = "/gpfs/group/dfc13/default/dasha/mlarson/pm_trial/nutau/*.i3.bz2"
# for filename in glob.glob(data_file_1):
#     file_list_n.append(filename)
# for filename in glob.glob(data_file_2):
#     file_list_n.append(filename)
# for filename in glob.glob(data_file_3):
#     file_list_n.append(filename)

data_file_1 = "/gpfs/group/dfc13/default/dasha/mlarson/pm_trial/corsika/*.i3.bz2"
for filename in glob.glob(data_file_1):
    file_list_c.append(filename)

# data_file_1 = "/gpfs/group/dfc13/default/dasha/mlarson/L2/data/*.i3.bz2"
# for filename in glob.glob(data_file_1):
#     file_list_d.append(filename)

#Create a ST configuration service that will be used by the icetray module.

i_frame = geofile.pop_frame()
g_frame = geofile.pop_frame()
geometry = g_frame["I3Geometry"]
polygon = []
top = 0
bottom = 0
count = 0
for i, j in geometry.omgeo:
    if (i.om == 39) and ((i.string == 26) or (i.string == 27) or (i.string == 37) or (i.string == 46) or (i.string == 45) or (i.string == 35)):
        polygon.append([i.string, j.position.x, j.position.y])

    if (i.om == 60) and (i.string == 36):
        bottom = j.position.z

    if (i.om == 39) and (i.string == 36):
        top = j.position.z

path = mplPath.Path(np.array([[polygon[0][1], polygon[0][2]], [polygon[1][1], polygon[1][2]], [polygon[2][1], polygon[2][2]], [polygon[3][1], polygon[3][2]], [polygon[4][1], polygon[4][2]], [polygon[5][1], polygon[5][2]]]))

now=datetime.now()
print now

values = []
none = []

VertexFits = ['MPEFit', 'SplineMPE_MPEFit_0_0', 'SplineMPE_MPEFit_Vertex_1_1', 'SplineMPE_MPEFit_VertexSRT_1_1', 'SplineMPE_MPEFit_VertexFid_1_1','SplineMPE_MPEFit_VertexSRTFid_1_1', 'SplineMPE_Primary_0_0']

#PmNames = ['SplitInIcePulses_BestFineFit_prob_obs_0s_cherdat_1','SplitInIcePulses_PrimaryOrig_prob_obs_0s_cherdat_1', 'SplitInIcePulses_SplineMPE_Primary_0_0_prob_obs_0s_cherdat_1']

PmNames0 = ['SplitInIcePulses_PrimaryOrig_prob_obs_0s_cherdat_1']

PmNames1 = ['SplitInIcePulses_BestFineFits_0_prob_obs_0s_cherdat_1','SplitInIcePulses_BestFineFits_1_prob_obs_0s_cherdat_1','SplitInIcePulses_BestFineFits_2_prob_obs_0s_cherdat_1','SplitInIcePulses_BestFineFits_3_prob_obs_0s_cherdat_1','SplitInIcePulses_BestFineFits_4_prob_obs_0s_cherdat_1'] 

def Genie(frame,name):
    if name == "genie":
        return True
    else:
        return False
 
def event_counter(frame):                                                                      
    global count                                                                                 
#        if "LineFit" not in frame.keys():                                                           
#            print count #frame["LineFit"]                                                          
    count += 1                                                                                
    if count%1000 ==0:
        print count
    

def Values(frame, year,  PmNames0, PmNames1):
    global values 
    vt = {}
    cor = {}
    arr = []

    if frame.Has("Weight") and frame.Has("VT_Q"):    
            
#        arr.append(frame["Energy"].value)
#        arr.append(frame["Contained"].value)
        arr.append(frame["Weight"].value)
        arr.append(frame["Prim_Q"].value)
        arr.append(frame["Prim_WQ"].value)
        arr.append(frame["Prim_P"].value)
        arr.append(frame["Prim_H"].value)
        arr.append(frame["Prim_Tot"].value)
        arr.append(np.array(frame["VT_Q"]))
        arr.append(np.array(frame["VT_WQ"]))
        arr.append(np.array(frame["VT_P"]))
        arr.append(np.array(frame["VT_H"]))
        arr.append(np.array(frame["VT_Tot"]))
        arr.append(np.array(frame["CT_Q"]))
        arr.append(np.array(frame["CT_WQ"]))
        arr.append(np.array(frame["CT_P"]))
        arr.append(np.array(frame["CT_H"]))
        arr.append(np.array(frame["CT_Tot"]))
        
        values.append(arr)
        if len(values)%100 == 0:
            print len(values)/100
  
        return True

    print"Something Missing"

def MaxCorTrack(frame):
    cor = set()
    Qs = []
    WQs = []
    Ps = []
    Hs = []
    Tot = []
    name_prim_p = "TrackHits_PrimaryOrig_SplitInIcePulses_coincObsProbsList_1"
    name_prim_q = "TrackHits_PrimaryOrig_SplitInIcePulses_coincObsQsList_1"
    prim_p = frame[name_prim_p]
    prim_q = frame[name_prim_q]
    p_oms = []
    for om, pr in prim_p:
        if pr:
            p_oms.append(om)

    prim_coinc =[]

    for k in frame.keys():
        qs = []
        wqs = []
        lqs = []
        hs = 0
        primc = 0
        if ("TrackHits_CorridorTrack" in k) and ("coincObsQsList" in k):
            trk_n = k.split("_")[2]
            if not trk_n in cor:
                cor.update([trk_n])
                namep = "TrackHits_CorridorTrack_{0}_SplitInIcePulses_coincObsProbsList_1".format(trk_n)
                nameq = "TrackHits_CorridorTrack_{0}_SplitInIcePulses_coincObsQsList_1".format(trk_n)
                probs =  frame[namep]
                chs =  frame[nameq]
                tot_doms = len(chs)

                #non_zero = {}
                for om, q in chs:
                    if len(q) != 0:
                        if om in p_oms:
                            primc = primc +1
                        
                        qs.append(sum(q))
                        hs = hs + 1
                        lqs.append(len(q))
                        wqs.append(np.array(probs[om]).dot(np.array(q)))

                if hs != 0:
                    Qs.append(sum(qs))
                    WQs.append(sum(wqs)) 
                    Ps.append(sum(lqs))
                    Hs.append(hs)
                    Tot.append(tot_doms)
                    prim_coinc.append(primc)
                    

    frame["CT_Q"] = dataclasses.I3VectorDouble(Qs)
    frame["CT_WQ"] = dataclasses.I3VectorDouble(WQs)
    frame["CT_P"] = dataclasses.I3VectorDouble(Ps)
    frame["CT_H"] = dataclasses.I3VectorDouble(Hs)
    frame["CT_Tot"] = dataclasses.I3VectorDouble(Tot)
    frame["CT_PC"] = dataclasses.I3VectorDouble(prim_coinc)

def MaxVetoTrack(frame):
    cor = set()
    Qs = []
    WQs = []
    Ps = []
    Hs = []
    Tot = []
    name_prim_p = "TrackHits_PrimaryOrig_SplitInIcePulses_coincObsProbsList_1"
    name_prim_q = "TrackHits_PrimaryOrig_SplitInIcePulses_coincObsQsList_1"
    prim_p = frame[name_prim_p]
    prim_q = frame[name_prim_q]
    pqs = []
    pwqs = []
    plqs = []
    phs = 0
    p_oms = []
    for om, pr in prim_p:
        if pr:
            p_oms.append(om)
            pqs.append(sum(prim_q[om]))
            phs = phs + 1
            plqs.append(len(prim_q[om]))
            pwqs.append(np.array(pr).dot(np.array(prim_q[om])))

    frame["Prim_Q"] = dataclasses.I3Double(sum(pqs))
    frame["Prim_WQ"] = dataclasses.I3Double(sum(pwqs))
    frame["Prim_P"] = dataclasses.I3Double(sum(plqs))
    frame["Prim_H"] = dataclasses.I3Double(phs)
    frame["Prim_Tot"] = dataclasses.I3Double(len(prim_p))


    prim_coinc = []

    for k in frame.keys():
        qs = []
        wqs = []
        lqs = []
        hs = 0
        primc = 0
        if ("TrackHits_VetoFit" in k) and ("coincObsQsList" in k):
            trk_n = k.split("_")[2]
            if not trk_n in cor:
#                print trk_n
                namep = "TrackHits_VetoFit_{0}_SplitInIcePulses_coincObsProbsList_1".format(trk_n)
                nameq = "TrackHits_VetoFit_{0}_SplitInIcePulses_coincObsQsList_1".format(trk_n)
                probs =  frame[namep]
                chs =  frame[nameq]
                cor.update([trk_n])
                tot_doms = len(chs)

                for om, pr in probs:
                    if pr:
                        if om in p_oms:
                            primc = primc +1
                        
                        qs.append(sum(chs[om]))
                        hs = hs + 1
                        lqs.append(len(chs[om]))
                        wqs.append(np.array(pr).dot(np.array(chs[om])))
                if qs:
                    prim_coinc.append(primc)
                    Qs.append(sum(qs))
                    WQs.append(sum(wqs)) 
                    Ps.append(sum(lqs))
                    Hs.append(hs)
                    Tot.append(tot_doms)

#    print prim_coinc
    frame["VT_Q"] = dataclasses.I3VectorDouble(Qs)
    frame["VT_WQ"] = dataclasses.I3VectorDouble(WQs)
    frame["VT_P"] = dataclasses.I3VectorDouble(Ps)
    frame["VT_H"] = dataclasses.I3VectorDouble(Hs)
    frame["VT_Tot"] = dataclasses.I3VectorDouble(Tot)
    frame["VT_PC"] = dataclasses.I3VectorDouble(prim_coinc)

def Energy(frame):

    vertex = 0
    energy = 0
    time = 0
    direc = 0
    contained = False
    if frame.Has("I3MCTree"):
        prim = dataclasses.get_most_energetic_neutrino(frame['I3MCTree'])
        vertex = prim.pos
        direc = prim.dir
        time = prim.time
        energy = prim.energy
    else:
        print "No Primary"
        return
            
    if vertex > 0 and energy > 0:   
        if (vertex.z < top) and (vertex.z > bottom) and path.contains_point((vertex.x, vertex.y)):
            contained = True
    
    frame["Energy"] = dataclasses.I3Double(energy)
    frame["Time"] = dataclasses.I3Double(time)
    frame["Dir"] = dataclasses.I3Direction(direc)
    frame["Vertex"] = dataclasses.I3Position(vertex)
    frame["Contained"] = icetray.I3Bool(contained)

    return

#@icetray.traysegment
def TestCuts(name, file_list, year):
  
    tray = I3Tray()
    tray.AddModule("I3Reader", name +"reader", FilenameList = file_list)
    
#    tray.AddModule(event_counter, 'counter', Streams = [icetray.I3Frame.Physics])        
    #tray.AddModule(Energy, name +"GE", If = lambda f: Genie(f,name))   
    tray.AddModule(MaxCorTrack, name +"CT")   
    tray.AddModule(MaxVetoTrack, name +"VT")   
    tray.Add(Values, name + "getvalues", year = year, PmNames0 = PmNames0, PmNames1 = PmNames1)
    #tray.AddModule('I3Writer', 'writer', Filename='TestCutsNeutrino2'+'.i3.bz2', Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics], DropOrphanStreams=[icetray.I3Frame.DAQ])
    tray.AddModule('TrashCan','thecan')
    tray.Execute()
    tray.Finish()

    return

print count
import pickle
TestCuts("genie", file_list = file_list_c, year = "12")
values_c = values[:]
name_f = "Corsika_TH_prim.pkl"
#name_f = "Genie_TH_numu_{0}_{1}.pkl".format(p1,p2)
#np.save(name_f,values_n)
output = open(name_f,"wb")
pickle.dump(values_c, output, -1)
output.close()
#none_n = none[:]
del values[:]
#del none[:]
now=datetime.now()
print "corsika finished", now, len(values_c)

