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

p = sys.argv[1]
#p2 = sys.argv[2]
#print p, int(p)   
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

data_file_1 = "/gpfs/group/dfc13/default/dasha/mlarson/output_1/nue/Pm_genie_ic.12640.0000{:02d}*.i3.bz2".format(int(p))
#data_file_1 = "/gpfs/group/dfc13/default/dasha/mlarson/output_1/numu/Pm_genie_ic.14640.0000{:02d}*.i3.bz2".format(int(p))
#print data_file_1
#data_file_1 = "/gpfs/group/dfc13/default/dasha/mlarson/output_1/nutau/Pm_genie_ic.16640.0000{:02d}*.i3.bz2".format(int(p))
# data_file_1 = "/gpfs/group/dfc13/default/dasha/mlarson/output_1/nue/*.i3.bz2"
# data_file_2 = "/gpfs/group/dfc13/default/dasha/mlarson/output_1/numu/*.i3.bz2"
# data_file_3 = "/gpfs/group/dfc13/default/dasha/mlarson/output_1/nutau/*.i3.bz2"
for filename in glob.glob(data_file_1):
    file_list_n.append(filename)
print file_list_n    
# for filename in glob.glob(data_file_2):
#     file_list_n.append(filename)
# for filename in glob.glob(data_file_3):
#     file_list_n.append(filename)

# data_file_1 = "/gpfs/group/dfc13/default/dasha/mlarson/output_1/corsika/*.i3.bz2"
# for filename in glob.glob(data_file_1):
#     file_list_c.append(filename)

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
    
def TrackStats(frame):    
    corrs = set()
    vfits = set()
    for k in frame.keys():  
        if ("CorridorTrack" in k):
            words = k.split("_")
            idx = words.index("CorridorTrack")+1
            name = "CorridorTrack"+"_"+words[idx]
            if name not in corrs:
                corrs.add(name)
  
        if ("VetoFit" in k):
            words = k.split("_")
            idx = words.index("VetoFit")+1
            name = "VetoFit"+"_"+words[idx]
            if name not in vfits:
                vfits.add(name)

    tracknames_c =  list(corrs)
    tracknames_v = list(vfits)
    trk_arr = []
    corr_arr = []
    for trackname in tracknames_v:
        logl = 0
        ilogl = 0
        flogl = 0
        tlogl = 0
        pm = 0
        prbs = []
        nprbs = []
        qs = []
        nqs = []

        track = frame[trackname]
        track_arr = [track.pos.x,track.pos.y,track.pos.z,track.dir.zenith,track.dir.azimuth]
        fitname = "SplineMPE_"+trackname+"_1_1"
        fit_arr = []
        if frame.Has(fitname):
            fit = frame[fitname]
            fit_arr = [fit.pos.x,fit.pos.y,fit.pos.z,fit.dir.zenith,fit.dir.azimuth]

        for k in frame.keys():
            if ("SplineMPE" in k) and (trackname in k) and ("FitParams" in k):
                tlogl = frame[k].logl
            if ("LLHCalcMPE" in k) and (trackname in k):
                logl = frame[k].logl
            if ("LLHCalcInitMPE" in k) and (trackname in k):
                ilogl = frame[k].logl
            if ("LLHCalcFinMPE" in k) and (trackname in k):
                flogl = frame[k].logl    
            if ("prob_obs_0s" in k) and (trackname in k):
                pm = frame[k].value
            if ("coincObsQsList" in k) and (trackname in k):
                Ps = frame[k]
                for om, value in Ps:
                    if value:
                        qs.append(value[0])
            if ("noncoincObsQsList" in k) and (trackname in k):
                Ps = frame[k]
                for om, value in Ps:
                    if value:
                        nqs.append(value[0])
            
                        
            if ("coincObsProbsList" in k) and (trackname in k):
                Probs = frame[k]
                for om, value in Probs:
                    if value:
                        prbs.append(value[0])

            if ("noncoincObsProbsList" in k) and (trackname in k):
                Probs = frame[k]
                for om, value in Probs:
                    if value:
                       nprbs.append(value[0])
            

        trk_arr.append([logl,pm,ilogl,tlogl,flogl,prbs,qs,nprbs,nqs,track_arr,fit_arr])
        
    # for trackname in tracknames_c:
    #     logl = 0
    #     pm = 0
    #     prbs = []
    #     qs = []

    #     track = frame[trackname]
    #     track_arr = [track.pos.x,track.pos.y,track.pos.z,track.dir.zenith,track.dir.azimuth]
       
    #     for k in frame.keys():
            
    #         if ("LLHCalcMPE" in k) and (trackname in k):
    #             logl = frame[k].logl
    #         if ("prob_obs_0s" in k) and (trackname in k):
    #             pm = frame[k].value
    #         if ("coincObsQsList" in k) and (trackname in k):
    #             Ps = frame[k]
    #             for om, value in Ps:
    #                 if value:
    #                     qs.append(value[0])
    #         if ("coincObsProbsList" in k) and (trackname in k):
    #             Probs = frame[k]
    #             for om, value in Probs:
    #                 if value:
    #                     prbs.append(value[0])

    #     corr_arr.append([logl,pm,prbs,qs,track_arr])

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
    
    arr = [energy, contained, time, vertex.x, vertex.y, vertex.z, direc.azimuth, direc.zenith]
    
    global values
    #print "SSSSS", arr+[trk_arr,corr_arr]
#    values.append(arr+[trk_arr,corr_arr])
    values.append(arr+trk_arr)
    if len(values)%100 == 0:
        print len(values)/100


#@icetray.traysegment
def TestCuts(name, file_list, year):
  
    tray = I3Tray()
    tray.AddModule("I3Reader", name +"reader", FilenameList = file_list)
    
    tray.AddModule(TrackStats, name +"CT")   
    #tray.AddModule('I3Writer', 'writer', Filename='TestCutsNeutrino2'+'.i3.bz2', Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics], DropOrphanStreams=[icetray.I3Frame.DAQ])
    tray.AddModule('TrashCan','thecan')
    tray.Execute()
    tray.Finish()

    return

#print count
import pickle
TestCuts("genie", file_list = file_list_n, year = "13")
values_n = values[:]
#name_f = "Genie_FullTest.pkl"
name_f = "Genie_TH2_nue_{0}.npy".format(p)
#np.save(name_f,values_n)
output = open(name_f,"wb")
pickle.dump(values_n, output, -1)
output.close()
# #none_n = none[:]
# del values[:]
# #del none[:]
# now=datetime.now()
print "genie finished", now, len(values_n)

