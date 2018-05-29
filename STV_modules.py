import numpy
import pickle
from icecube import phys_services,linefit,gulliver,gulliver_modules,spline_reco,StartingTrackVeto
from I3Tray import I3Units
from icecube import dataclasses, icetray
from icecube import dataio, phys_services
from I3Tray import *
import copy


def RemoveFromFrame(frame,names):
    for name in names:
        if frame.Has(name):
            del frame[name]

def GetLLH(frame,FitName):
    geometry = frame["I3Geometry"]
    log_calc = spline_reco.I3SplineRecoLikelihood()
    log_calc.SetGeometry(geometry)
    log_calc.SetEvent(frame)
    logl = log_calc.GetLogLikelihood(I3EventHypothesis(frame['FitName']))
    name = 'SplineRecoLLH_{0:s}'.format("FitName")


def DoGetLLH(tray, name, FitNames):
    for fitname in FitNames:
        tray.AddSegment(GetLLH,"GetLLH"+fitname,
                        FitName=fitname)
    
    
def DoSplineReco(tray,name,Pulses,Seed,Llh,LogName,Spline,AngStep,DistStep,If=lambda frame: True):
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

def NSegmentVector(frame,FitName,N=1):
    if N%2==0:
        print "n=",N,"is even! Change this!"
        sys.exit(910)
    try:
        basep=copy.deepcopy(frame[FitName])
    except:
        return True

    #shift to closest approach to 0,0,0    
    origin_cap=phys_services.I3Calculator.closest_approach_position(basep,dataclasses.I3Position(0,0,0))
    
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
                for om in veto_hits:
                    newfit = copy.deepcopy(frame[CoGFit])                                             
                    newfit.dir=dataclasses.I3Direction(newfit.pos-geo[om[0]].position)
                    name = 'VetoFit_{0:02d}{1:02d}'.format(om[0].string,om[0].om)
                    frame[name]=newfit
                return True

        if frame.Has("CoGTime") and frame.Has(SafeFit):
            fit = copy.deepcopy(frame[SafeFit])                                                
            if fit.fit_status == dataclasses.I3Particle.OK:  
                for om in veto_hits:
                    newfit = copy.deepcopy(frame[SafeFit])                                           
                    newfit.pos=frame["CoGPos"]
                    newfit.dir=dataclasses.I3Direction(newfit.pos-geo[om[0]].position)
                    name = 'VetoFit_{0:02d}{1:02d}'.format(om[0].string,om[0].om)
                    frame[name]=newfit
                return True

#    print "Can't make VetoFits"
    return False

def MakeCorridorFits(frame, CoGFit, SafeFit):
    input_file = open('CorridorTracks.pkl', 'rb')
    data = pickle.load(input_file)
    if len(data) != 0:
        
        if frame.Has(CoGFit):
            fit = copy.deepcopy(frame[CoGFit])                                                
            if fit.fit_status == dataclasses.I3Particle.OK:  
                for i in data:
                    newfit = copy.deepcopy(frame[CoGFit])                                             
                    newfit.dir=dataclasses.I3Direction(newfit.pos-dataclasses.I3Position(i[1]))
                    name = i[0]
                    frame[name]=newfit
                return True

        if frame.Has("CoGTime") and frame.Has(SafeFit):
            fit = copy.deepcopy(frame[SafeFit])                                                
            if fit.fit_status == dataclasses.I3Particle.OK:
                for i in data:
                    newfit = copy.deepcopy(frame[SafeFit])                                            
                    newfit.dir=dataclasses.I3Direction(newfit.pos-dataclasses.I3Position(i[1]))
                    name = i[0]
                    frame[name]=newfit
                return True

    print "Can't make CorridorFits"
    return False

    
def CascadeDir(frame, CoGFit, SRTPulsesFid):                                                   
    
    geo = frame["I3Geometry"].omgeo                                                             
    if frame.Has(SRTPulsesFid) and frame.Has(CoGFit):
        fit=copy.deepcopy(dataclasses.I3Particle(frame[CoGFit]))
        if not fit.fit_status == dataclasses.I3Particle.OK:  
            return
        pulsesf = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,SRTPulsesFid)               
        if  len(pulsesf) == 0:                                                                      
            return                                                                                  
        max_time = 0                                                                                
        max_dist = 0                                                                               
        max_dist_om = 0
        max_time_om = 0

        for om, pulseSeries in pulsesf:                                                             
            for pulse in pulseSeries:                                                               
                dist = (fit.pos-geo[om].position).magnitude                               
                time = pulse.time - fit.time
                if dist > max_dist:
                    max_dist = dist
                    max_dist_om = om
                if time > max_time:
                    max_time = time
                    max_time_om = om
        
        if (max_dist_om != 0):
            fit_dist =  copy.deepcopy(dataclasses.I3Particle(frame[CoGFit]))
            fit_dist.dir=dataclasses.I3Direction(geo[max_dist_om].position-fit_dist.pos)
            frame["FarthestDistFit"]=fit_dist
            frame['FarthestDistFitDOM']=dataclasses.I3VectorOMKey([max_dist_om])
        
        if (max_time_om != 0):
            fit_time =  copy.deepcopy(dataclasses.I3Particle(frame[CoGFit]))
            fit_time.dir=dataclasses.I3Direction(geo[max_time_om].position-fit_time.pos)
            frame['LatestTimeFit']=fit_time
            frame['LatestTimeFitDOM']=dataclasses.I3VectorOMKey([max_time_om])
               
    return           


def CoGMedIC(frame, PulsesFid):                                                                      
    geometry = frame["I3Geometry"]                                                                   
    if frame.Has(PulsesFid):                                                                         
        pulsesf = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,PulsesFid)                       
        if  len(pulsesf) == 0:                                                                       
            return                                                                                   
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
    return           

def GetBestCoarseFit(frame, FitNames):

    input_file = open('CorridorTracks.pkl', 'rb')
    data = pickle.load(input_file)
    names = [i[0] for i in data]
    poss = [i[1] for i in data]

    geo = frame["I3Geometry"].omgeo
    rlogl = []
    rlogl_best = 9999
    fitname_best = " "
    best_fits = []

    for fitname in FitNames:
        if frame.Has(fitname):
            for k in frame.keys():
                if ("FitParams" in k) and (fitname in k):
                    rlogl.append([frame[k].rlogl,fitname])

    if len(rlogl) == 0:
        print "Coarse Fits failed, can't get the best"
        return False
    
    rlogl_sort = sorted(rlogl, key = lambda rlogl: rlogl[0])

    for fit in rlogl_sort:
        if (fit[0]-rlogl_sort[0][0])/rlogl_sort[0][0]*100 < 2:
            best_fits.append(fit[1])
    
    if len(best_fits) < 1:
        print "Fit probably failed"
        return False

    if len(best_fits) > 1:
        best_fits = best_fits[:1] 

    for i in range(0,len(best_fits)):
        
        coarse_fit = copy.deepcopy(frame[best_fits[i]])
        name_parts = best_fits[i].split("_")
        if name_parts[1] == 'VetoFit':
            om = OMKey(int(name_parts[2][:2]),int(name_parts[2][2:4]))
            pos =geo[om].position
        elif name_parts[1] == 'FarthestDistFit':
            om = frame['FarthestDistFitDOM'][0]
            pos =geo[om].position
        elif name_parts[1] == 'CorridorTrack':
            idx = names.index('CorridorTrack_'+name_parts[2])
            pos = poss[idx]
        elif name_parts[1] == 'LatestTimeFit':
            om = frame['LatestTimeFitDOM'][0]
            pos =geo[om].position
        else:
            pos = coarse_fit.pos-(coarse_fit.dir*700)
        
        frame['BestCoarseFits_'+str(i)] = coarse_fit
        frame['BestCoarseFits_'+str(i)+"FitParams"] = copy.deepcopy(frame[best_fits[i]+"FitParams"])
        frame['BestCoarseFits_'+str(i)+'Pos'] = dataclasses.I3VectorI3Position([coarse_fit.pos, pos])
        

def MakeFineFitPos(frame, CoarseFits):

    geo = frame["I3Geometry"].omgeo
    for fit in CoarseFits:
        if frame.Has(fit):
        
            N0= frame[fit+'Pos'][0]
            N1= frame[fit+'Pos'][1]
            
            P0 = numpy.array([N0.x,N0.y,N0.z])
            P1 = numpy.array([N1.x,N1.y,N1.z])

            U = P1-P0 #direction of the coarse track                                                           
            norm = numpy.linalg.norm(U)
            U = U / norm

            V =  numpy.array([0,-U[2],U[1]])
            norm = numpy.linalg.norm(V)
            V = V / norm

            W =  numpy.array([U[1]**2+U[2]**2, -U[0]*U[1], -U[0]*U[2]])
            norm = numpy.linalg.norm(W)
            W = W / norm
  
            field_fid = [[-20,21],[-20,21]]
            field_veto = [[-60,61],[-60,61]]
            step_fid = [10,10]
            step_veto = [30,30]
        
            V_fid = numpy.arange(field_fid[0][0],field_fid[0][1],step_fid[0])
            W_fid = numpy.arange(field_fid[1][0],field_fid[1][1],step_fid[1])
            V_veto = numpy.arange(field_veto[0][0],field_veto[0][1],step_veto[0])
            W_veto = numpy.arange(field_veto[1][0],field_veto[1][1],step_veto[1])

            VW_fid = [[v,w] for v in V_fid for w in W_fid]
            VW_veto = [[v,w] for v in V_veto for w in W_veto]

            pos_fid = []
            pos_veto = []

            for pos in VW_fid:
                shift = pos[0]*V+pos[1]*W
                pos_fid.append(P0+shift)
        
            for pos in VW_veto:
                shift = pos[0]*V+pos[1]*W
                pos_veto.append(P1+shift)
                    
            save1 = []
            save2 = []

            for i in range(0,len(pos_fid)):
                save1.append(dataclasses.I3Position(pos_fid[i][0],pos_fid[i][1],pos_fid[i][2]))
            frame[fit+'GridFid'] = dataclasses.I3VectorI3Position(save1)


            for i in range(0,len(pos_veto)):
                save2.append(dataclasses.I3Position(pos_veto[i][0],pos_veto[i][1],pos_veto[i][2]))
            frame[fit+'GridVeto'] = dataclasses.I3VectorI3Position(save2)

def MakeFineFits(frame,CoarseFits):

    for coarse_fit in CoarseFits:
        if frame.Has(coarse_fit):
            pos_fid = frame[coarse_fit+'GridFid']
            pos_veto = frame[coarse_fit+'GridVeto']
    
            for i in range(0, len(pos_fid)):
                for j in range(0, len(pos_veto)):
                    fit= copy.deepcopy(dataclasses.I3Particle(frame[coarse_fit]))
                    fit.pos = pos_fid[i]
                    fit.dir=dataclasses.I3Direction(fit.pos-pos_veto[j])
                    index = i*len(pos_fid)+j
                    number = CoarseFits.index(coarse_fit)
                    name = 'FineFit_{0:02d}_{1:04d}'.format(number,index)
                    frame[name]=fit

def GetBestFineFit(frame, FitNames):

    rlogl = []
    rlogl_best = 9999
    fitname_best = " "
    best_fits =[] 

    for fitname in FitNames:
        if frame.Has(fitname):
            for k in frame.keys():
                if ("FitParams" in k) and (fitname in k):
                    rlogl.append([frame[k].rlogl,fitname])
                    
    if len(rlogl) == 0:
        print "Fine Fits failed, can't get the best"
        return False
    
    rlogl_sort = sorted(rlogl, key = lambda rlogl: rlogl[0])
    fine_fit_name = rlogl_sort[0][1]
    
    for fit in rlogl_sort:
        if (fit[0]-rlogl_sort[0][0])/rlogl_sort[0][0]*100 < 2:
            best_fits.append(fit[1])
    
    if len(best_fits) < 1:
        print "Fit probably failed"
        return False

    if len(best_fits) > 5:
        best_fits = best_fits[:5] 

    for i in range(0,len(best_fits)):
        frame['BestFineFits_'+str(i)] = copy.deepcopy(frame[best_fits[i]])
        frame['BestFineFits_'+str(i)+"FitParams"] = copy.deepcopy(frame[best_fits[i]+"FitParams"])

#    print best_fits
    fine_fit = copy.deepcopy(frame[fine_fit_name])
    fine_fit_params = copy.deepcopy(frame[fine_fit_name+"FitParams"])

    frame['BestFineFit'] = fine_fit
    frame['BestFineFitNames'] = dataclasses.I3VectorString(best_fits)
    frame['BestFineFitFitParams'] = fine_fit_params
            
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
    
