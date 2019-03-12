##D.Pankova 02/28/2019
##If you are running on L2, you might want to use precuts to remove the obvious muons first

from __future__ import division
import icecube
from icecube import icetray, dataclasses, dataio, phys_services, STTools
from I3Tray import *
import copy

#Save truth information
def Primary(frame,FitName,GenieOrCorsika="Genie"):
    if frame.Has(FitName):
        fit = copy.deepcopy(frame[FitName])
        if not fit.fit_status == dataclasses.I3Particle.OK:
            print "Bad Primary Seed"
            return False

        if frame.Has("I3MCTree"):
            if GenieOrCorsika=="Genie":
                prim = dataclasses.get_most_energetic_neutrino(frame['I3MCTree'])
            elif GenieOrCorsika=="Corsika":    
                prim = dataclasses.get_most_energetic_muon(frame['I3MCTree'])
            else:
                print "Specify if neutrino or muons"
                return False
            fit.pos = prim.pos
            fit.dir = prim.dir                
            fit.time = prim.time
            frame['Primary'] = fit
            return True

    print "No Primary"
    return False

#Get 4 variables, which I found have the strongest effect on background removal
#Number of STR cleaned pulses in DeepCore, Number of pulses in Veto
#Veto to Deep core charge ratio and Vertex Z position estimation
def CalculateVars(frame, SRTPulses, SRTPulsesFid, SRTPulsesVeto, PulsesVeto):                               
    geometry = frame["I3Geometry"]
    vertex_z = -999
    first_hit_t = 1e6
    h_fid = 0                                                                           
    p_veto = 0     
    c_fid = 0                                                                          
    c_veto = 0                                                                              
 
    if frame.Has(SRTPulsesFid): 
        for omkey in frame[SRTPulsesFid]:
            h_fid += 1
            for pulse in omkey[1]:
                c_fid += pulse.charge
        if (h_fid == 0) or (c_fid == 0):           
            print "CalclateVars: nothing SRT in DeepCore"
            return False
    else:
        print "CalclateVars: no SRTPulsesFid"
        return False

    if frame.Has(SRTPulsesVeto): 
        for omkey in frame[SRTPulsesVeto]:
            for pulse in omkey[1]:
                c_veto += pulse.charge   
    else:
        print "CalclateVars: no SRTPulsesVeto"
        return False

    if frame.Has(PulsesVeto):
        for omkey in frame[PulsesVeto]:  
            for pulse in omkey[1]:
                p_veto += 1
    else:
        print "CalclateVars: no PulsesVeto"
        return False

    if frame.Has(SRTPulses):            
        for omkey in frame[SRTPulses].apply(frame):
            for pulse in omkey[1]:
                if pulse.time < first_hit_t:
                    first_hit_t = pulse.time
                    vertex_z = geometry.omgeo[omkey[0]].position.z
    else:
        print "CalclateVars: no Pulses"
        return False
               
    c_ratio = c_veto/c_fid 
     
    frame["Vars_Hits_Fid"] = dataclasses.I3Double(h_fid)                          
    frame["Vars_Pulses_Veto"] = dataclasses.I3Double(p_veto)                         
    frame["Vars_Charge_Ratio"] = dataclasses.I3Double(c_ratio)                      
    frame["Vars_Vertex_Z"] = dataclasses.I3Double(vertex_z)
    print "Vars:", h_fid, p_veto, c_ratio, vertex_z

    #Functions below take into account correlations between some variables
    func_1 = p_veto - (0.3*h_fid**2 + 40)/(-h_fid) - 50
    frame["Vars_Func_1"] = dataclasses.I3Double(func_1)
    func_2 = c_ratio + 0.0007*(func_1 + 34)*(func_1 - 6)
    frame["Vars_Func_2"] = dataclasses.I3Double(func_2)
    print "Func: ", func_1, func_2

    return True
#Do the cut
def PreCut(frame):
    if frame.Has("Vars_Func_1") and frame.Has("Vars_Charge_Ratio") and frame.Has("Vars_Vertex_Z"):
        func_1 = frame["Vars_Func_1"].value
        cr = frame["Vars_Charge_Ratio"].value
        vz = frame["Vars_Vertex_Z"].value
        if (func_1 < 0) and (cr < 0.2) and (vz < -190):
            return True
        else:
            print "Didn't pass PreCut"
            return False
    else:
        print "Didn't have something in PreCut"
        return False

#Remove events with a lot of compatible hits
def THCut(frame, Pulses, FitNames, Threshold):
    for fitname in FitNames:
        for k in frame.keys():
            if ("TrackHits_{0}_{1}".format(fitname,Pulses) in k) and ("coincObsQsList" in k):
                lists = frame[k]
                Qsum = []
                for om, val in lists: #Find non zero hits
                    Qsum.append(sum(val))
                SQ = sum(Qsum)     
                if SQ> Threshold:
                    print "Didn't Pass TrackHits Cut"
                    print "Q ={0}, Threshold = {1}".format(SQ,Threshold)
                    return False
    return True           
