##D.Pankova 02/28/2019 IceCube
##Functions used in Starting Track Veto for Low energies

import numpy
import pickle
from icecube import phys_services, linefit, gulliver, gulliver_modules, spline_reco, StartingTrackVetoLE, TrackHits
from I3Tray import I3Units
from icecube import dataclasses, icetray
from icecube import dataio, phys_services
from I3Tray import *
from operator import itemgetter
import copy
    
def DoSplineReco(tray,name,Pulses,Seed,Llh,LogName,Spline,AngStep,DistStep,
                 If=lambda frame: True):
    
    # setup minuit, parametrization, and bayesian priors for general use                            
    tray.AddService("I3GulliverMinuitFactory", "Minuit%s" % (Llh) + LogName,
                    Algorithm="SIMPLEX",
                    MaxIterations=100000,
                    Tolerance=0.1,
                    )

    tray.AddService("I3SimpleParametrizationFactory", "SimpleTrack%s" % (Llh) + LogName,
                    StepX = DistStep*I3Units.m,
                    StepY = DistStep*I3Units.m,
                    StepZ = DistStep*I3Units.m,
                    StepT = 3.333333  *I3Units.ns,
                    StepZenith = AngStep*I3Units.degree,
                    StepAzimuth= AngStep*I3Units.degree,
                    BoundsX = [-2000*I3Units.m, 2000*I3Units.m],
                    BoundsY = [-2000*I3Units.m, 2000*I3Units.m],
                    BoundsZ = [-2000*I3Units.m, 2000*I3Units.m],
                    )

    tray.AddService("I3BasicSeedServiceFactory", "SplineSeed%s" % (Llh) + LogName,
                    FirstGuesses=[Seed],
                    )

    tray.AddService("I3SplineRecoLikelihoodFactory","LLHSpline%s" % (Llh)+LogName,
                    PhotonicsService=Spline,
                    Pulses=Pulses,
                    Likelihood=Llh,
                    #NoiseRate=10*I3Units.hertz,                                                     
                    )

    tray.AddModule( "I3SimpleFitter", "Spline%s" % (Llh) + LogName,
                    OutputName = "Spline%s" % (Llh) + LogName,
                    SeedService="SplineSeed%s" % (Llh) + LogName,
                    Parametrization="SimpleTrack%s" % (Llh) + LogName,
                    LogLikelihood="LLHSpline%s" % (Llh) + LogName,
                    Minimizer="Minuit%s" % (Llh) + LogName,
                    If= lambda frame: frame.Has(Seed),
                    )

def EvalLLH(tray, name, Llh, Pulses, Spline, FitNames):
    #Evaluate Likelihood given a track and pulses
    for fitname in FitNames:  
        tray.AddService("I3SplineRecoLikelihoodFactory","LLHSplineEval%s"%(Llh)+"_"+fitname,
                        PhotonicsService=Spline,
                        Pulses=Pulses,
                        Likelihood=Llh,
                        #NoiseRate=10*I3Units.hertz,          
                        )

        tray.AddModule("I3LogLikelihoodCalculator", "LLHCalc%s" % (Llh)+"_" + fitname,
                        FitName=fitname,
                        LogLikelihoodService="LLHSplineEval%s" % (Llh)+"_"+fitname,
                        )


def DoVetoFits(tray, name, Pulses,Llh, FitNames, Spline, AngStep, DistStep):
    #Run a fitter on given seed tracks
    for fitname in FitNames:
        tray.AddSegment(DoSplineReco,"DoSplineReco"+fitname,
                        Pulses=Pulses,
                        Seed=fitname,
                        Llh=Llh,
                        LogName = "_" + fitname,
                        Spline=Spline,
                        AngStep=AngStep,
                        DistStep=DistStep)

def DoVetoPulseFits(tray, name, Llh, FitNames, Spline, AngStep, DistStep):
    #Run a fitter on given seed tracks and corrseponsing pulse series
    for fitname in FitNames:
        tray.AddSegment(DoSplineReco,"DoSplineReco"+fitname,
                        Pulses=fitname+"_Pulses",
                        Seed=fitname,
                        Llh=Llh,
                        LogName = "_" + fitname,
                        Spline=Spline,
                        AngStep=AngStep,
                        DistStep=DistStep)


def SelectLLH(frame, TrackNames, N):
    #Select most promising Tracks (the rest will be deleted)
    #N with minimum LLH and N with maximum number of compatible Hits
    lists = []
    fitnames = []
    for fitname in TrackNames:
        if frame.Has(fitname):
            logl = 10**10
            p = 0
            q = 0
            for k in frame.keys():
                if ("LLHCalc" in k) and (fitname in k):
                    if not numpy.isnan(frame[k].logl):
                        logl = frame[k].logl
                    
                if ("TrackHits" in k) and ("coincObsQsList" in k) and (fitname in k):
                    Ps = frame[k]
                    for om, charges in Ps:
                        if not len(charges) == 0:
                            p = p + 1
                            q = q + sum(charges)
            lists.append([logl,p,q,fitname])

    ps = sorted(lists, key=itemgetter(1), reverse=True)
    qs = sorted(lists, key=itemgetter(2), reverse=True)
    ps = ps[:N] 
    logls = sorted(lists, key=itemgetter(0))
    logls = logls[:N]
    frame["TrackHits_MaxCompHits"]=dataclasses.I3Double(ps[0][1])
    frame["TrackHits_MaxCompCharge"]=dataclasses.I3Double(qs[0][2])
    print "TrackHits_MaxCompHits = {0:.3f}".format(ps[0][1])
    print "TrackHits_MaxCompCharge = {0:.3f}".format(qs[0][2])
    
    for it in logls:
        fitnames.append(it[3])

    for it in ps:
        if not (it[3] in fitnames):
            fitnames.append(it[3])
            
    for fitname in TrackNames:
        if not (fitname in fitnames):
            for k in frame.keys():
                if (fitname in k):
                    del frame[k]
                
            

def NSegmentVector(frame,FitName,N=1):
    #Make Segments out of tracks, Required for STV (made by K.Jero)
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
        #print "Put", FitName+"_"+str(N)+"_segments", "in the frame"                             
        del segments

def DoSTV(tray, name, Pulses, FitNames, NSeg, Spline, MinCADist, DistType):
    #Run STV
    for fitname in FitNames:
        #Create vectors for STV                                                             
        tray.Add(NSegmentVector,"NSegmentVector_"+fitname+"_"+str(NSeg),
                 FitName=fitname,
                 N=NSeg)

        tray.Add("StartingTrackVetoLE","STV_"+fitname+"_"+str(NSeg),
                 Pulses=Pulses,
                 Photonics_Service=Spline,
                 Miss_Prob_Thresh=1.1, #No cuttting on Pm here
                 Fit=fitname,
                 Particle_Segments=fitname+"_"+str(NSeg)+"_segments",
                 Distance_Along_Track_Type=DistType,
                 Supress_Stochastics=False,
                 Cascade = True,
                 Norm = False,
                 Min_CAD_Dist=MinCADist)

def DoTrackHits(tray, name, Pulses, FitNames, NSeg, Spline, MinCADist):
    #Run TrackHits Module, looking for compatible Hits
    for fitname in FitNames:
        tray.Add(NSegmentVector,"NSegmentVectorTH_"+fitname+"_"+str(NSeg),
                 FitName=fitname,
                 N=NSeg)        
        
        tray.Add("TrackHits","TH_"+fitname+"_"+str(NSeg),
                 Pulses=Pulses,
                 Photonics_Service=Spline,
                 Percent=0.01, #Same as inside STV for LE, no point 
                 Fit=fitname,  #of having a different one here
                 Particle_Segments=fitname+"_"+str(NSeg)+"_segments",
                 Min_CAD_Dist=MinCADist)
        
      

def MakeVetoTracks(frame, CoGFit, SafeFit, PulsesVeto, CoGPosName, CoGTimeName):
    #connect one hit in the Veto to a Vertex, CoG fit vertex or CoG in DeepCore
    geo = frame["I3Geometry"].omgeo
    if frame.Has(PulsesVeto):
        veto_hits = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,PulsesVeto)
        if  len(veto_hits) == 0:
            print "No Hits in Veto!"
            return False

        #Use CoGFit Vertex, if CoGfit exists
        if frame.Has(CoGFit):
            fit = copy.deepcopy(frame[CoGFit])                                         
            if fit.fit_status == dataclasses.I3Particle.OK:
                num = 0
                for om in veto_hits:
                    newfit = copy.deepcopy(frame[CoGFit])                              
                    newfit.dir=dataclasses.I3Direction(newfit.pos-geo[om[0]].position)
                    name = 'VetoFit_{0:05d}'.format(num)
                    frame[name]=newfit
                    num = num +1
                return True

        #If not take CoG as a vertex    
        if frame.Has(CoGTimeName) and frame.Has(SafeFit):
            fit = copy.deepcopy(frame[SafeFit])                                            
            if fit.fit_status == dataclasses.I3Particle.OK:
                num = 0
                for om in veto_hits:
                    newfit = copy.deepcopy(frame[SafeFit])                                  
                    newfit.pos=frame[CoGPosName]
                    newfit.time=frame[CoGTimeName].value
                    newfit.dir=dataclasses.I3Direction(newfit.pos-geo[om[0]].position)
                    name = 'VetoFit_{0:05d}'.format(num)
                    frame[name]=newfit
                    num = num +1
                return True

    print "Can't make VetoFits"
    return False

    
def MakeVetoPulses(frame, PulsesVeto, PulsesFid):
    #Make pulse series consisting of fid pulses and one hit in veto
    cor = set() #to make sure there are no repeats
    geo = frame["I3Geometry"].omgeo
    if not frame.Has(PulsesFid) or not frame.Has(PulsesVeto):        
        print "Can't make VetoPulses"
        return False
    
    veto_hits = copy.deepcopy(dataclasses.I3RecoPulseSeriesMap.from_frame(frame,PulsesVeto))
    if len(veto_hits) == 0:
        print "MakeVetoPulses: No Hits in Veto!"
        return False
    
    #look for TrackHits output and read out compatiable DOMS
    for k in frame.keys():
        if ("TrackHits_VetoFit" in k) and ("coincObsQsList" in k):
            trk = k.split("_")[2]
            if trk in cor:
                continue

            Qs = copy.deepcopy(frame[k])
            cor.update([trk])
            veto_OMs = []

            #Find non zero hits
            for om, q in Qs: 
                if sum(q) != 0:
                    veto_OMs.append(om)

            total_hits = copy.deepcopy(dataclasses.I3RecoPulseSeriesMap.from_frame(frame,PulsesFid))
            if len(total_hits) == 0:
                print "MakeVetoPulses: No Hits in DC!"
                return False

            #Go through all compatiable hits
            if veto_OMs:
                for om in veto_OMs: 
                    if om in total_hits:
                        print "MakeVetoPulses: Veto and Fid pulses not separated!"
                        return False
                    #Make Pulses Series out of them
                    total_hits[om] = veto_hits[om]
     
            name = 'VetoFit_{0}_Pulses'.format(trk)
            frame[name]=total_hits
                        
    return True

#Calculate and book Pmiss
def Pmiss(frame, Pulses, TrackNames):
    LLHs = []
    Pms = []
    Qs = []
    for fitname in TrackNames:
        for k in frame.keys():
            if ("TrackHits_{0}_{1}".format(fitname,Pulses) in k) and ("coincObsQsList" in k):
                lists = frame[k]
                Qsum = []
                for om, val in lists: #Find non zero hits
                    Qsum.append(sum(val))
                Qs.append(sum(Qsum))
            if ("LLHCalc" in k) and (fitname in k):   
                LLHs.append(frame[k].logl)
    
            if ("prob_obs_0s" in k) and (fitname in k):
                Pms.append(frame[k].value)
    
    zipped = zip(LLHs,Qs,Pms)
    maxQ = numpy.max(Qs)
    maxQ_LLH_list = [i[0] for i in zipped if i[1] == maxQ]
    maxQ_minLLH = numpy.min(maxQ_LLH_list)
    pm_best = [i[2] for i in zipped if i[1] == maxQ and i[0] == maxQ_minLLH]
    pm_max = numpy.max(Pms)
    pm_min = numpy.min(Pms)
    pm_mean = numpy.mean(Pms)

    frame["Pm_Max_STV"] = dataclasses.I3Double(pm_max)
    frame["Pm_Min_STV"] = dataclasses.I3Double(pm_min)
    frame["Pm_Mean_STV"] = dataclasses.I3Double(pm_mean)
    frame["Pm_Best_STV"] = dataclasses.I3Double(pm_best[0])
    print "Pm_Max = {0:.3E}, Pm_Min = {1:.3E}".format(pm_max, pm_min)
    print "Pm_Mean = {0:.3E}, Pm_Best = {1:.3E}".format(pm_mean, pm_best[0])
    return True

#Find Center of Gravity of the Deep Core pulses (approximate vertex)
def CoGMedIC(frame, PulsesFid, CoGPosName, CoGTimeName):                            
    geometry = frame["I3Geometry"]                                          
    if frame.Has(PulsesFid):                                                      
        pulsesf = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,PulsesFid)    

        if  len(pulsesf) == 0:     
            print "CoGMedIC: Pulses are empty"
            return False                                                               

        cog_time = []                                                      
        cog_x = []                                                    
        cog_y = []                     
        cog_z = []                                                         

        for om, pulseSeries in pulsesf:                                  
            for pulse in pulseSeries:                                            
                cog_x.append(geometry.omgeo[om].position.x)                        
                cog_y.append(geometry.omgeo[om].position.y)                             
                cog_z.append(geometry.omgeo[om].position.z)                         

        cog_pos  = dataclasses.I3Position(numpy.median(cog_x), numpy.median(cog_y), numpy.median(cog_z))       
        cog_time = []                                                                 
        distance = 0                                                     
         
        for om, pulseSeries in pulsesf:                                           
            distance = (cog_pos - geometry.omgeo[om].position).magnitude                
            for pulse in pulseSeries:                                         
                cor_time = abs(pulse.time - distance/dataclasses.I3Constants.c_ice)      
                cog_time.append(cor_time)

        cog_t = numpy.median(cog_time)                                             
        frame[CoGTimeName] = dataclasses.I3Double(cog_t)                                    
        frame[CoGPosName] = dataclasses.I3Position(cog_pos)                               

    else:
        print "CoGMedIC: No Fid Pulses"
        return False          

#Functions below are not used in the current version od STV for LE
#But may still be useful at some point
def Intersection(cog, mid, p1, p2):
    #Find if two lines intesect for making corridor tracks
    #First Line: Connect CoG of Event in Deep Core and a position(x,y) 
    #in the middle of the outer edge of the corridor (outer edge of IC)
    #Second Line: Take positions(x,y) of two strings at the 
    #inner edge of the coridor near Deep Core

    tr = mid-cog
    cor = p2-p1
    denom = tr.x*cor.y-cor.x*tr.y
    #never colliner or parallel
    #if denom == 0 : return None # collinear
    denom_is_positive = denom > 0
    top = cog-p1
    s_numer = tr.x*top.y-tr.y*top.x
    if (s_numer < 0) == denom_is_positive:
        return False # no collision
    t_numer = cor.x*top.y-cor.y*top.x
    if (t_numer < 0) == denom_is_positive:
        return False # no collision
    if (s_numer > denom) == denom_is_positive or (t_numer > denom) == denom_is_positive:
        return False # no collision
    #If collision is detected
    #Find the intersection point and make sure it's not too close
    #To the edge of the Second Line, i.e. a muon can pass through
    #the corridor without getting too close to any strings
    t = t_numer/denom
    intsec = dataclasses.I3Position(cog.x+(t*tr.x), cog.y+(t*tr.y),0)
    d1 = (intsec-p1).magnitude
    d2 = (intsec-p2).magnitude
    if (d1 < 30) or (d2 < 30):
        return False
    return True

#Put tracks through the "active" corridors
def MakeCorridorFits(frame, CoGFit, SafeFit):
    input_file = open('Corridors.pkl', 'rb')
    data = pickle.load(input_file)
    if len(data) == 0:
        print "MakeCorridorFits: No File"
        return False

    #find CoG
    fit = 0
    if frame.Has(CoGFit):
        fit = copy.deepcopy(frame[CoGFit])
        if fit.fit_status == dataclasses.I3Particle.OK:
            #make tracks        
            for cor in data:
                use_cor = Intersection(fit.pos, data[cor][0][0], 
                                       data[cor][0][1], data[cor][0][2]) 
                if use_cor == True:
                    for trk in data[cor][1]:
                        newfit = copy.deepcopy(fit)                                  
                        newfit.dir = dataclasses.I3Direction(
                            newfit.pos - trk[1])
                        frame[trk[0]]=newfit
            return True
        
   
    if frame.Has(SafeFit) and frame.Has("CoGPos"):   
        fit = copy.deepcopy(frame[SafeFit])
        if fit.fit_status == dataclasses.I3Particle.OK:
            fit.pos = frame["CoGPos"]
            fit.time = frame["CoGTime"].value
            #make tracks        
            for cor in data:
                use_cor = Intersection(fit.pos, data[cor][0][0], 
                                       data[cor][0][1], data[cor][0][2]) 
                if use_cor == True:
                    for trk in data[cor][1]:
                        newfit = copy.deepcopy(fit)                                  
                        newfit.dir = dataclasses.I3Direction(
                            newfit.pos - trk[1])
                        frame[trk[0]]=newfit
            return True
    
    print "MakeCorridorFits: No good fits"
    return False


        
            
    
