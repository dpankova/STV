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

data_file_1 = "/gpfs/group/dfc13/default/dasha/mlarson/pm_trial/nue/*.i3.bz2"
data_file_2 = "/gpfs/group/dfc13/default/dasha/mlarson/pm_trial/numu/*.i3.bz2"
data_file_3 = "/gpfs/group/dfc13/default/dasha/mlarson/pm_trial/nutau/*.i3.bz2"
for filename in glob.glob(data_file_1):
    file_list_n.append(filename)
for filename in glob.glob(data_file_2):
    file_list_n.append(filename)
for filename in glob.glob(data_file_3):
    file_list_n.append(filename)

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

PmNames0 = ['SplitInIcePulses_PrimaryOrig_prob_obs_0s_cherdat_1', 'SplitInIcePulses_SplineMPE_Primary_0_0_prob_obs_0s_cherdat_1']

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
    poss = []
    pms0 = []
    pms1 = []
    arr = []

    if frame.Has("Energy") and frame.Has("Contained") and frame.Has("Weight"):    
            
        for pm in PmNames0:
            if not frame.Has(pm):
                return False
            pms0.append(frame[pm].value)

        print pms0
        for pm in PmNames1:
            if frame.Has(pm):
                pms1.append(frame[pm].value)
        
        print pms1
        arr.append(pms0)
        arr.append(pms1)
        arr.append(frame["Energy"].value)
        arr.append(frame["Contained"].value)
        arr.append(frame["Weight"].value)
        values.append(arr)

        done_fits = False
        for k in frame.keys():
            if "SplineMPE_VetoFit" in k:
                done_fits = True

        arr.append(done_fits)
        if len(values)%1000 == 0:
            print len(values)/1000
  
        return True

    print"Something Missing"

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
    tray.AddModule(Energy, name +"GE", If = lambda f: Genie(f,name))   
    tray.Add(Values, name + "getvalues", year = year, PmNames0 = PmNames0, PmNames1 = PmNames1)
    #tray.AddModule('I3Writer', 'writer', Filename='TestCutsNeutrino2'+'.i3.bz2', Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics], DropOrphanStreams=[icetray.I3Frame.DAQ])
    tray.AddModule('TrashCan','thecan')
    tray.Execute()
    tray.Finish()

    return

print count
TestCuts("genie", file_list = file_list_n, year = "13")
values_n = values[:]
np.save("Genie_Pmiss_M1_Corr_Mult.npy",np.array(values_n))
#none_n = none[:]
del values[:]
#del none[:]
now=datetime.now()
print "genie finished", now, len(values_n)

# TestCuts("corsika", file_list = file_list_n, year = "12")
# values_n = values[:]
# none_n = none[:]
# np.save("/storage/home/dup193/work/estes/new_estes/Corsika_BDTVals.npy",np.array(values_n))
# del values[:]
# del none[:]
# now=datetime.now()
# print "corsika1 finished", now, len(values_n)

# TestCuts("data", file_list = file_list_d, year = "13")
# values_d = values[:]
# none_d = none[:]
# np.save("Data_BDTVals_mySRT.npy",np.array(values_d))
# del values[:]
# del none[:]
# now=datetime.now()
# print "data finished", now, len(values_d)

#print 'total (N C D) =', len(values_n), len(values_c), len(values_d) 
#print 'no pulses (NCD) =', len(none_n), len(none_c), len(none_d) 
#pickle.dump(solar_n, open("genie_solar_MedIC.pkl","wb"), protocol = pickle.HIGHEST_PROTOCOL)
