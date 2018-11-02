import numpy
import pickle
from icecube import phys_services, linefit, gulliver, gulliver_modules, spline_reco, StartingTrackVeto, TrackHits
from I3Tray import I3Units
from icecube import dataclasses, icetray
from icecube import dataio, phys_services
from I3Tray import *
from operator import itemgetter
import copy


def RemoveFromFrame(frame,names):
    for name in names:
        if frame.Has(name):
            del frame[name]
    
    
def DoSplineReco(tray,name,
                 Pulses,Seed,
                 Llh,LogName,
                 Spline,AngStep,DistStep,
                 If=lambda frame: True):
    # setup minuit, parametrization, and bayesian priors for general use                            
    tray.AddService("I3GulliverMinuitFactory", "Minuit%s" % (Llh) + LogName,
                    Algorithm="SIMPLEX",
                    MaxIterations=100000,
                    Tolerance=100000000000000000,
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
    for fitname in FitNames:
#        tray.AddModule(lambda frame: frame.Has(fitname),"CheckEval%s"%(Llh)+"_"+fitname)  
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


def DoVetoFits(tray, name, Pulses,Llh,LogName,FitNames, Spline, AngStep, DistStep):
    for fitname in FitNames:
        tray.AddSegment(DoSplineReco,"DoSplineReco"+fitname,
                        Pulses=Pulses,
                        Seed=fitname,
                        Llh=Llh,
                        LogName = "_" + fitname +LogName,
                        Spline=Spline,
                        AngStep=AngStep,
                        DistStep=DistStep)

def DoVetoPulseFits(tray, name, Llh,LogName,FitNames, Spline, AngStep, DistStep):
    for fitname in FitNames:
        tray.AddSegment(DoSplineReco,"DoSplineReco"+fitname,
                        Pulses=fitname+"_Pulses",
                        Seed=fitname,
                        Llh=Llh,
                        LogName = "_" + fitname +LogName,
                        Spline=Spline,
                        AngStep=AngStep,
                        DistStep=DistStep)


def SelectLLH(frame, TrackNames, N):
    lists = []
    fitnames = []
    for fitname in TrackNames:
        if frame.Has(fitname):
            logl = 10**10
            p = 0
            for k in frame.keys():
                if ("LLHCalc" in k) and (fitname in k):
                    if not numpy.isnan(frame[k].logl):
                        logl = frame[k].logl
                    
                if ("coincObsQsList" in k) and (fitname in k):
                    Ps = frame[k]
                    for om, value in Ps:
                        if value:
                            p = p + 1

            lists.append([logl,p,fitname])
    #print lists        
    ps = sorted(lists, key=itemgetter(1), reverse=True)
    ps = ps[:N] 
    logls = sorted(lists, key=itemgetter(0))
    logls = logls[:N]
    #print "p", ps 
    #print "l", logls
    for it in logls:
        fitnames.append(it[2])

    for it in ps:
        if not (it[2] in fitnames):
            fitnames.append(it[2])
    #print "f", fitnames        
    for fitname in TrackNames:
        if not (fitname in fitnames):
            for k in frame.keys():
                if (fitname in k):
                    del frame[k]
                
            

def NSegmentVector(frame,FitName,N=1):
    if frame.Has(FitName):

        if N%2==0:
            print "n=",N,"is even! Change this!"
            sys.exit(910)
        try:
            basep=copy.deepcopy(frame[FitName])
        except:
            return True

    #shift to closest approach to 0,0,0    
        origin_cap = phys_services.I3Calculator.closest_approach_position(
            basep,dataclasses.I3Position(0,0,0))
    
        basep_shift_d=numpy.sign(origin_cap.z - basep.pos.z) *\
            numpy.sign(basep.dir.z) *\
            (origin_cap-basep.pos).magnitude
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

        RemoveFromFrame(frame,[FitName+"_"+str(N)+"_segments"])
        frame[FitName+"_"+str(N)+"_segments"]=dataclasses.I3VectorI3Particle(segments)
    #    print "Put", FitName+"_"+str(N)+"_segments", "in the frame"                             
        del segments

def PrintPm(frame, FitNames):
    pm = 0
    rlogl = 0
    for fitname in FitNames:
        if frame.Has(fitname):
            for k in frame.keys():
                if ("prob_obs_0s" in k) and (fitname in k):
                    pm = frame[k].value
                if ("FitParams" in k) and (fitname in k):
                    rlogl = frame[k].rlogl
            print "{0} Rlogl = {1:.3e} Pm = {2:.3e}".format(fitname, rlogl, pm)
            print "Getting out", fitname, frame[fitname]


def DoSTV(tray, name, Pulses, FitNames, NSeg, PmCut, Spline, MinCADist, DistType):
    for fitname in FitNames:
        
        #Create vectors for STV                                                             
        tray.Add(NSegmentVector,"NSegmentVector_"+fitname+"_"+str(NSeg),
                 FitName=fitname,
                 N=NSeg)
        
        #STV                                                                                
        tray.Add("StartingTrackVeto","STV_"+fitname+"_"+str(NSeg),
                 Pulses=Pulses,
                 Photonics_Service=Spline,
                 Miss_Prob_Thresh=PmCut,
                 Fit=fitname,
                 Particle_Segments=fitname+"_"+str(NSeg)+"_segments",
                 Distance_Along_Track_Type=DistType,
                 Supress_Stochastics=False,
                 Cascade = True,
                 Norm = False,
                 Min_CAD_Dist=MinCADist)

def PrintCoincQP(frame, Pulses, FitName, NSeg):
    ps = 0
    qs = 0
    name_p = "TrackHits_"+FitName+"_"+Pulses+"_"+"coincObsPsList_"+str(NSeg)
    name_q = "TrackHits_"+FitName+"_"+Pulses+"_"+"coincObsQsList_"+str(NSeg)
    if frame.Has(name_p):
        for k in frame.keys():
            if (name_p in k):
                ps = sum(frame[k])
            if (name_q in k):
                qs = sum(frame[k])    
        print "{0} Qs = {1:.3e} Ps = {2:.3e}".format(FitName, qs, ps)


def DoTrackHits(tray, name, Pulses, FitNames, NSeg, Percent, Spline, MinCADist):
    for fitname in FitNames:
        #Create vectors for STV  
        #print "NSegmentVectorTH_",fitname,str(NSeg)
        tray.Add(NSegmentVector,"NSegmentVectorTH_"+fitname+"_"+str(NSeg),
                 FitName=fitname,
                 N=NSeg)        
        #STV                                                                                
        tray.Add("TrackHits","TH_"+fitname+"_"+str(NSeg),
                 Pulses=Pulses,
                 Photonics_Service=Spline,
                 Percent=Percent,
                 Fit=fitname,
                 Particle_Segments=fitname+"_"+str(NSeg)+"_segments",
                 Min_CAD_Dist=MinCADist)
        
        # tray.Add(PrintCoincQP,"THPrint_"+fitname+"_"+str(NSeg),
        #          Pulses = Pulses,
        #          FitName=fitname,
        #          NSeg = NSeg)
      

def MakeVetoFits(frame, CoGFit, SafeFit, PulsesVeto):
    geo = frame["I3Geometry"].omgeo
    if frame.Has(PulsesVeto):
        veto_hits = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,PulsesVeto)
        if  len(veto_hits) == 0:
            print "No Hits in Veto!"
            return False

        if frame.Has(CoGFit):
            fit = copy.deepcopy(frame[CoGFit])                                                
            if fit.fit_status == dataclasses.I3Particle.OK:
                num = 0
                for om in veto_hits:
                    newfit = copy.deepcopy(frame[CoGFit])                                             
                    newfit.dir=dataclasses.I3Direction(newfit.pos-geo[om[0]].position)
#                    name = 'VetoFit_{0:02d}{1:02d}'.format(om[0].string,om[0].om)
                    name = 'VetoFit_{0:05d}'.format(num)
                    frame[name]=newfit
                    num = num +1
                return True

        if frame.Has("CoGTime") and frame.Has(SafeFit):
            fit = copy.deepcopy(frame[SafeFit])                                                
            if fit.fit_status == dataclasses.I3Particle.OK:
                num = 0
                for om in veto_hits:
                    newfit = copy.deepcopy(frame[SafeFit])                                           
                    newfit.pos=frame["CoGPos"]
                    newfit.time=frame["CoGTime"].value
                    newfit.dir=dataclasses.I3Direction(newfit.pos-geo[om[0]].position)
                   # name = 'VetoFit_{0:02d}{1:02d}'.format(om[0].string,om[0].om)
                    name = 'VetoFit_{0:05d}'.format(num)
                    frame[name]=newfit
                    num = num +1
                return True

#    print "Can't make VetoFits"
    return False

    
def MakeVetoPulses(frame, PulsesVeto, PulsesFid):
    #Make pulse series = fid pulses +one hit in veto
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
#    print veto_hits
    for k in frame.keys():
        if ("TrackHits_VetoFit" in k) and ("coincObsQsList" in k):

            trk = k.split("_")[2]
            if trk in cor:
                continue

            Qs = copy.deepcopy(frame[k])
            cor.update([trk])
            veto_OMs = []
            
            for om, q in Qs: #Find non zero hits
                if sum(q) != 0:
                    veto_OMs.append(om)

            if not veto_OMs: #Didn't find any non-zero hits
                continue
#            print veto_OMs

            total_hits = copy.deepcopy(dataclasses.I3RecoPulseSeriesMap.from_frame(frame,PulsesFid))
#            print "total", total_hits

            if len(total_hits) == 0:
                print "MakeVetoPulses: No Hits in DC!"
                return False
    
            for om in veto_OMs: #Go through all compatiable hits
                if om in total_hits:
                    print "MakeVetoPulses: Veto and Fid pulses not separated!"
#                    print "OM", om
                    return False
                #Make Pulses Series out of them
                total_hits[om] = veto_hits[om]

            name = 'VetoFit_{0}_Pulses'.format(trk)
            frame[name]=total_hits
                        
    return True
    

def Intersection(cog, mid, p1, p2):

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

    # collision detected
    t = t_numer/denom

    intsec = dataclasses.I3Position(cog.x+(t*tr.x), cog.y+(t*tr.y),0)
    d1 = (intsec-p1).magnitude
    d2 = (intsec-p2).magnitude
    
    if (d1 < 30) or (d2 < 30):
        return False
    return True


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

def CoGMedIC(frame, PulsesFid):                                                    
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
        frame["CoGTime"] = dataclasses.I3Double(cog_t)                                     
        frame["CoGPos"] = dataclasses.I3Position(cog_pos)                               
    else:
        print "CoGMedIC: No Fid Pulses"
        return False          

def CheckTH(frame, NHitsCut):
    cor = set()
    Qs = [] #Charge
    Ps = [] #Pulses
    Hs = [] #Hits
    Tot = [] #Total #DOMs around the track
    MaxNHits = 0 #Max number of hits in one track, 
                 #the variable to cut on

    for k in frame.keys():
        qs = []
        wqs = []
        lqs = []
        hs = 0
        #Get TH data from the frame
        if ("TrackHits_VetoFit" in k) and ("coincObsQsList" in k):
            trk_n = k.split("_")[2]
            if not trk_n in cor:
                namep = "TrackHits_VetoFit_{0}_SplitInIcePulses_coincObsProbsList_1".format(trk_n)
                nameq = "TrackHits_VetoFit_{0}_SplitInIcePulses_coincObsQsList_1".format(trk_n)
                probs =  frame[namep]
                chs =  frame[nameq]
                cor.update([trk_n])
                tot_doms = len(chs)

                #Summ charge pulses for each om, if nonzero
                for om, q in cks:
                    if q:
                        qs.append(sum(q))
                        hs = hs + 1
                        lqs.append(len(q))

                #Sum charge for each track    
                if qs:
                    Qs.append(sum(qs))
                    Ps.append(sum(lqs))
                    Hs.append(hs)
                    Tot.append(tot_doms)

    if not Hs:
        MaxNHits = 0
    else:
        Hs.sort(reverse=True)
        MaxNHits = Hs[0]

    
        
            
def RLogLCut(frame,fits_to_try=[]):
    for f in fits_to_try:
        if frame.Has(f):
            fps=frame[f+"FitParams"]
            rlogl=fps.rlogl
            print f,rlogl
            if numpy.isnan(rlogl):
                return False
        else:
            return False
    return True
    
