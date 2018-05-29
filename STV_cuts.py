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

def CountAndCRatio(frame, PulsesFid, PulsesVeto):                                             
    if frame.Has(PulsesFid) and frame.Has(PulsesVeto):                                                
        c_fid_count = 0                                                                               
        c_vet_count = 0                                                                               
        p_count = 0                                                                                   
        for omkey in frame[PulsesFid]:                                                                
            for pulse in omkey[1]:                                                                    
                c_fid_count += pulse.charge                                                           
                p_count += 1                                                                          
                                                                                                      
        frame["NumPulsesFid"] = dataclasses.I3Double(p_count)                                      

        p_count = 0                                                                                   
        for omkey in frame[PulsesVeto]:                                                               
            for pulse in omkey[1]:                                                                    
                c_vet_count += pulse.charge                                                           
                p_count += 1                                                                          
                                                                                                      
        frame["NumPulsesVeto"] = dataclasses.I3Double(p_count)                         

        # if c_fid_count != 0:                                                                       
        #     frame["CRatio"+name] = dataclasses.I3Double(c_vet_count/c_fid_count)                   
        # else:                                                                                       
        #     frame["CRatio"+name] = dataclasses.I3Double(99)                                        
    return                    


def PulsesCut(frame,NumPFid,NumPVeto):
    if frame.Has("NumPulsesVeto") and frame.Has("NumPulsesFid"):
        fid = frame["NumPulsesFid"].value
        veto = frame["NumPulsesVeto"].value
        if (fid>NumPFid) and (veto<NumPVeto):
            return True

    print "Didn't pass Pulses Cut"
    return False

