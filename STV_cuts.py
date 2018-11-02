#!/usr/bin/python                                                                                   
from __future__ import division
import icecube
from icecube import icetray, dataclasses, dataio, phys_services, STTools
from I3Tray import *
import copy

def Primary(frame,FitName):

    # if frame.Has("I3MCTree"):
    #     I3MCT=frame["I3MCTree"]
    #     for prim in I3MCT.primaries:
    #         if (abs(prim.pdg_encoding)==12) or (abs(prim.pdg_encoding)==14) or (abs(prim.pdg_encoding)==16):
    #             fit = copy.deepcopy(frame[FitName])
    #             fit.pos = prim.pos
    #             fit.dir = prim.dir                
    #             fit.time = prim.time
    #             frame['Primary'] = fit
    #             frame['PrimaryOrig'] = fit
    #             return True
    if frame.Has(FitName):
        fit = copy.deepcopy(frame[FitName])
        if not fit.fit_status == dataclasses.I3Particle.OK:
            print "Bad Primary Seed"
            return False

        if frame.Has("I3MCTree"):
            prim = dataclasses.get_most_energetic_neutrino(frame['I3MCTree'])
            fit.pos = prim.pos
            fit.dir = prim.dir                
            fit.time = prim.time
            frame['Primary'] = fit
            frame['PrimaryOrig'] = fit
            return True
  

    print "No Primary"
    return False

def CalculateVars(frame, Pulses, PulsesFid, PulsesVeto):                                             
    if frame.Has(PulsesFid) and frame.Has(PulsesVeto) and frame.Has(Pulses):                                                
        geometry = frame["I3Geometry"]
        vertex_z = -999
        first_hit_t = 1e6
        h_fid = 0                                                                               
        p_veto = 0     
        c_fid = 0                                                                               
        c_veto = 0                                                                              
                                                                            
        for omkey in frame[PulsesFid]:
            h_fid += 1
            for pulse in omkey[1]:                                                                    
                c_fid += pulse.charge                                                           
        
        if h_fid == 0:
            print "CalclateVars: nothing in DeepCore"
            return False

        for omkey in frame[PulsesVeto]:                                                               
            for pulse in omkey[1]:                                                                    
                p_veto += pulse.charge  

        for omkey in frame[Pulses].apply(frame):
            for pulse in omkey[1]:
                if pulse.time < first_hit_t:
                    first_hit_t = pulse.time
                    vertex_z = geometry.omgeo[omkey[0]].position.z
        c_ratio = c_veto/c_fid            
        frame["Vars_Hits_Fid"] = dataclasses.I3Double(h_fid)                                                                        
        frame["Vars_Pulses_Veto"] = dataclasses.I3Double(p_veto)                         
        frame["Vars_Charge_Ratio"] = dataclasses.I3Double(c_ratio)                               
        frame["Vars_Vertex_Z"] = dataclasses.I3Double(vertex_z)
        print "Vars:", h_fid, p_veto, c_ratio, vertex_z
        
        func_1 = p_veto - (0.3*h_fid**2 + 40)/(-h_fid) - 50
        frame["Vars_Func_1"] = dataclasses.I3Double(func_1)

        func_2 = c_ratio + 0.0007*(func_1 + 34)*(func_1 - 6)
        frame["Vars_Func_2"] = dataclasses.I3Double(func_2)
        print "Func: ", func_1, func_2
    return

def PreCut(frame):
    if frame.Has("Vars_Func_1") and frame.Has("Vars_Func_2") and frame.Has("Vars_Vertex_Z"):
        func_1 = frame["Vars_Func_1"].value
        func_2 = frame["Vars_Func_2"].value
        vz = frame["Vars_Vertex_Z"].value
        if (func_1 < 0) and (func_2 < 1) and (vz < -190):
            return True
        else:
            print "Didn't pass PreCut"
            return False

