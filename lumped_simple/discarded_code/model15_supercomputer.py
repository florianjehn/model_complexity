#!/usr/bin/env python
# coding: utf-8 

"""

"""
from __future__ import division
import datetime

import cmf
import spotpy
from spotpy.parameter import Uniform as param
import sys
import os
import numpy as np
from datawriter_multi_objective import DataWriter


# Calibration time span
# 1979: Spin-Up
begin = 1980
end = 1985
# 1986 .. 1988: Validation

prefix='model15' ###### <------------- enter model name her 

# Number of runs
runs = 10


fnQ = 'GrebenauQTagMittel__1979_1990.txt'
fnT = 'Temp_max_min_avg_1979_1988.txt'
fnP = 'Prec_Grebenau_1979_1988.txt'


class Fulda_lumped(object):
    """
    Contains the whole model
    """
    def __init__(self,begin,end, with_valid_data = False, shift_one_day = False):
        """
        Initializes the model and builds the core setup  
        begin: begin of the calibration period
        eng: end of the calibration period
        with_calib_data: save also the data from the validation period
                        the calibration is still only done form 'begin' to 'end'
        """     
               # tr_S = Residence time of the water in the soil to the GW
        self.params = [param('tr_soil_GW',10.,400.),
                       # tr_soil_river = residence time from soil to river
                       param("tr_soil_fulda", 1.,170.),
                        # tr_surf = Residence time from surface 
                   #    param('tr_surf',0.001,30),
                       # tr_GW_l = residence time in the lower groundwate
                  #     param('tr_GW_l',1.,1000.),
                       # tr_GW_u = Residence time in the upper groundwater to the river
                       param('tr_GW_u_fulda',40.,350.),
                       # tr_GW_u_GW_l = residencete time to GW_l from GW_u
                #       param("tr_GW_u_GW_l", 10., 750.),
                       # tr_fulda = Residence time in the river (in days)
                 #      param('tr_fulda', 0., 3.5),  

                       # V0_soil = field capacity for the soil
                       param('V0_soil',5.,160.),

                        # beta_P_soil = Exponent that changes the form of the flux from the soil
                       param('beta_soil_GW',0.6,2.5),

                       # beta_fulda = exponent that changes the form of the flux from the soil 
                       param("beta_fulda", 2.0,5.5),

                       # ETV1 = the Volume that defines that the evaporation is lowered because of not enough water
                       param('ETV1',5.,100.),
                       # fETV0 = factor the ET is multiplied with when water is low
                       param('fETV0',0.,0.25),

                        # Rate of snow melt
                       param('meltrate',0.15,13.),
                       #  Snow_melt_temp = Temperature at which the snow melts (needed because of averaged temp
                       param('snow_melt_temp',0.0,3.5) ,
    
                       #Qd_max = maximal flux from lower groundwater to drinking water production
                       param('Qd_max', 0.,2.5),
                       # tw_thresholt = amount of water that can't be slurped out by the water pumps
                       param("TW_threshold", 50.,250.),

                       # LAI = leaf area index
                       param('LAI', 1.,12.),
                       # Canopy Closure
                       param("CanopyClosure",0.,0.5),

                       # Ksat = saturated conductivity of the soil 
                    #   param("Ksat", 0., 1)
                       ]        
        
        # loads the data  
        P,T,Tmin,Tmax,Q = self.loadPETQ()
        self.Q=Q
        # only use one core (quicker for small models)
        cmf.set_parallel_threads(1)
        # Generate a project with on ecell for a lumped model
        self.project = cmf.project()
        p = self.project
        
        # Add cell for soil and so on (x,y,z,area)
        c = p.NewCell(0,0,0,1000)
        
        # Add snow storage
        c.add_storage('Snow','S')
        cmf.Snowfall(c.snow,c)
        
        # Surfacewater is treated as a storage
       # c.surfacewater_as_storage()
        
        # Add the soil and groundwater layers to the soil cell
        soil = c.add_layer(2.0)
        gw_upper = c.add_layer(5.0) 
   #     gw_lower = c.add_layer(20.0)
        
        # Fill storages
        c.layers[0].volume = 15
        c.layers[1].volume = 80
      #  c.layers[2].volume = 120
    #    
        # Evapotranspiration
        cmf.HargreaveET(soil,c.transpiration)
        #cmf.PenmanMonteith()
        
        # Add the Fulda River
      #  self.fulda = p.NewOpenStorage(name="Fulda",x=0,y=0,z=0, area = 3.3*10**6)
        # Giving the Fulda a mean depth
     #   self.fulda.potential = 1.5   
   
        # add the drinking water outlet
        self.trinkwasser = p.NewOutlet('trinkwasser',20,0,0)
   
        # Outlet
        self.outlet = p.NewOutlet('outlet',10,0,0)
        
        # Storage for the interception
        I=c.add_storage('Canopy','C')
        
        # Rain
        self.makestations(P,T,Tmin,Tmax)
        self.project = p
        self.begin = begin
        self.end = end   
        self.with_valid_data = with_valid_data
        self.shift_one_day = shift_one_day
    
    def setparameters(self,
                      tr_soil_GW = 12.36870481, 
                      tr_soil_fulda = 12.,
                   #   tr_surf = 3.560855356,
                  #    tr_GW_l = 829.7188064, 
                      tr_GW_u_fulda = 270.05035, 
                  #    tr_GW_u_GW_l = 270., 
                   #   tr_fulda = 2.264612944,                     

                      V0_soil = 280.0850875,  
                      
                      beta_soil_GW=1.158865311, 
                      beta_fulda = 1.1,
                      
                      ETV1=2.575261852,
                      fETV0=0.014808919,
                      
                      meltrate = 4.464735097,
                      snow_melt_temp = 4.51938545,
                      
                      Qd_max = 0.250552812,
                      TW_threshold = 10.,
                      
                      LAI = 2.992013336,
                      CanopyClosure = 5.,
                      
                 #     Ksat = 0.02
                      ):  # this list has to be identical with the one above
        """
        sets the parameters, all parameterized connections will be created anew    
        """
        # Get all definitions from init method
        p = self.project
        c = p[0]
        outlet = self.outlet
    #    fulda = self.fulda
        trinkwasser = self.trinkwasser

        # Adjustment of the evapotranspiration
        c.set_uptakestress(cmf.VolumeStress(ETV1,ETV1 * fETV0))
        
        # Flux from the surfaces to the river
  #      cmf.kinematic_wave(c.surfacewater,fulda,tr_surf)
        # flux from surfaces to the soil (infiltration)
  #      cmf.SimpleInfiltration(c.layers[0], c.surfacewater) 

        # change the saturated conductivity of the soil
   #     c.layers[0].soil.Ksat = Ksat
         
        # Flux from soil to river (interflow)
        cmf.kinematic_wave(c.layers[0],outlet,tr_soil_fulda/V0_soil, V0 = V0_soil, exponent = beta_fulda)        
        # flux from the soil to the upper groundwater (percolation)
        cmf.kinematic_wave(c.layers[0], c.layers[1],tr_soil_GW, exponent=beta_soil_GW) 

        # flux from the upper groundwater to the river (baseflow)
        cmf.kinematic_wave(c.layers[1], outlet, tr_GW_u_fulda)               
        # flux from upper to lower groundwater (percolation)
 #       cmf.kinematic_wave(c.layers[1], c.layers[2],tr_GW_u_GW_l)#, exponent=beta_GW_u_GW_l) 
        
        # flux from the lower groundwater to river (baseflow)
    #    cmf.kinematic_wave(c.layers[2], fulda, tr_GW_l)        
        # Flux from the lower groundwater to the drinking water outlet
        # the fourths argument is the amount that is now allowed to be slurped 
        # out of the lower groundwater
        cmf.TechnicalFlux(c.layers[1],trinkwasser,Qd_max,TW_threshold,cmf.day)
        
        # Flux from drinking water to the river
        cmf.waterbalance_connection(trinkwasser, outlet)     
        
        # flux from the river to the outlet
    #    cmf.kinematic_wave(fulda, outlet, tr_fulda, exponent = beta_fulda) 
        
        # set snowmelt temperature
        cmf.Weather.set_snow_threshold(snow_melt_temp)        
        # Snowmelt at the surfaces
        snowmelt_surf = cmf.SimpleTindexSnowMelt(c.snow,c.surfacewater,c,rate=meltrate)

        # Splits the rainfall in interzeption and throughfall
        cmf.Rainfall(c.canopy,c, False, True)
        cmf.Rainfall(c.surfacewater,c, True, False)
        # Makes a overflow for the interception storage
        cmf.RutterInterception(c.canopy,c.surfacewater,c)
        # Transpiration on the plants is added
        cmf.CanopyStorageEvaporation(c.canopy,c.evaporation,c)
        # Sets the parameters for the interception       
        c.vegetation.LAI= LAI    
        # Defines how much throughfall there is (in %)
        c.vegetation.CanopyClosure = CanopyClosure
        
        
    def loadPETQ(self):
        """
        Loads climata and discharge data from the corresponding files fnQ, fnT and fnP              
        """    
        # Fixed model starting point
        begin = datetime.datetime(1979,1,1)
        step = datetime.timedelta(days=1)
        # empty time series
        P = cmf.timeseries(begin, step)
        P.extend(float(Pstr) for Pstr in open(fnP))
        
        Q = cmf.timeseries(begin,step)
        Q.extend(float(Qstr) for Qstr in open(fnQ))
        # Convert m3/s to mm/day
        Q *= 86400 * 1e3 / (2976.41 * 1e6)
        T = cmf.timeseries(begin,step)
        Tmin = cmf.timeseries(begin,step)
        Tmax = cmf.timeseries(begin,step)
        
        # Go through all lines in the file
        for line in open(fnT):
            columns = line.split('\t')
            if len(columns) == 3:
                Tmax.add(float(columns[0]))
                Tmin.add(float(columns[1]))
                T.add(float(columns[2]))
                
        return P,T,Tmin,Tmax,Q
        
    def makestations(self,P,T,Tmin,Tmax):
        """
        Creates the rainfall and the climate stations
        P = time series precipitation
        T, Tmin, Tmax = time series of mean temperatur, min and max         
        """
        rainstation = self.project.rainfall_stations.add('Grebenau avg',P,(0,0,0))
        self.project.use_nearest_rainfall()

        # Temperature data
        meteo = self.project.meteo_stations.add_station('Grebenau avg',(0,0,0))
        meteo.T = T
        meteo.Tmin = Tmin
        meteo.Tmax = Tmax
        self.project.use_nearest_meteo()
        
        return rainstation

    def runmodel(self,verbose=False):
        """
        starts the model
        if verboose = True --> give something out for every day    
        """
        try:
            # Creates a solver for the differential equations
            #solver = cmf.ImplicitEuler(self.project,1e-8)
            solver = cmf.CVodeIntegrator(self.project,1e-8)
            # usually the CVodeIntegrator computes the jakobi matrix only
            # partially to save computation time. But in models with low spatial
            # complexity this leads to a longer computational time
            # therefore the jakob matrix is computed completely to speed things up
            # this is done by LinearSolver = 0
            solver.LinearSolver = 0
            c = self.project[0]
            solver.max_step = cmf.h
            
            # New time series for model results (res - result)
            resQ = cmf.timeseries(self.begin,cmf.day)
            # starts the solver and calculates the daily time steps
            end = self.end
            if self.with_valid_data:
                end = datetime.datetime(1988,12,31)
            
            for t in solver.run(self.project.meteo_stations[0].T.begin, end,cmf.day):
                # Fill the results
                if t>=self.begin:
                    resQ.add(self.outlet.waterbalance(t))
                # Print a status report
                if verbose:
                    print (t,'Q=%5.3f, P=%5.3f' % (resQ[t],c.get_rainfall(t)))
                    
                # Print that one year was calculated, so one knows the model is still working
                #### comment this out if run on supercomputer to avoid spam ######
                #if t % cmf.year ==  cmf.year - cmf.year:
                 #   print("Finished one year")                    
                    
            # Return the filled result time series
            return resQ
        except RuntimeError:
            return np.array(self.Q[self.begin:self.end + datetime.timedelta(days=1)])*np.nan
            
    def simulation(self,vector):
        """
        SpotPy expects a method simulation. This methods calls setparameters
        and runmodels, so SpotPy is satisfied        
        """          
        
        paramdict = dict((pp.name,v) for pp,v in zip(self.params,vector))
        self.setparameters(**paramdict)
        resQ = self.runmodel()
        return np.array(resQ)


    def evaluation(self):
        """
        For Spotpy  
        """
        return np.array(self.Q[self.begin:self.end + datetime.timedelta(days=1)])
    

    def parameters(self):
        """
        For Spotpy  
        """ 
        return spotpy.parameter.generate(self.params)

    def objectivefunction(self,simulation,evaluation):
        """
        For Spotpy  
        """
        # to hit peaks better shift the timeseries by one day        
        if self.shift_one_day:
            simulation = simulation[:-1]
            evaluation = evaluation[1:]

        # if the validation data is added to the simulated data as well it should not
        # be used for calibration. To avoid this we have to shorten the list of 
        # the simulated data to the length of the calibration period            
        
        if self.with_valid_data:
            simulation = simulation[:len(evaluation)]
        logNS = spotpy.objectivefunctions.lognashsutcliff(evaluation, simulation)
        
        # calulate pbias here instead of problems to avoid problems with 
        # older spotpy versions
        sim = np.array(simulation)
        obs = np.array(evaluation)
        pbias = 100 * (float(np.sum( sim - obs )) / float(np.sum( obs )) )        
        
        rmse = spotpy.objectivefunctions.rmse(evaluation,simulation)
        standart_dev = obs.std()
        # rsr = Ratio between the root mean square error and the standart
        #    deviation of the measured data (see Moriasi et al 2007)
        rsr = rmse / standart_dev
        
#        print("logNS: "+str(logNS))
#        print("pbias: "+str(pbias))
#        print("rsr: "+str(rsr))
#        print()

        return [logNS, pbias, rsr]


if __name__ == '__main__': 
    # Import algorithm
    from spotpy.algorithms import lhs as Sampler

    # Find out if the model should run parallel (for supercomputer)
    parallel = 'mpi' if 'OMPI_COMM_WORLD_SIZE' in os.environ else 'seq'
    
    # Create the spotted model    
    model = Fulda_lumped(datetime.datetime(begin,1,1),
                          datetime.datetime(end,12,31), with_valid_data = True,
                        shift_one_day = True)
    if 'i' in sys.argv:
        runs = 0
    elif 'v' in sys.argv:
        sys.argv.remove('v')
        best = eval(open(prefix + '-best.dict').read())
        best.pop('Eff')
        model.setparameters(**best)
        model.begin = datetime.datetime(1986,1,1)
        model.end = datetime.datetime(1988,12,31)
        resQ = np.array(model.runmodel())
        model.plotvalidation(np.array(resQ))
        runs = 0
    elif len(sys.argv)>1:
        runs = int(sys.argv[1])
    if runs:
        sampler = Sampler(model, parallel=parallel)
      #  sampler.datawriter = DataWriter(prefix,model.params, model.begin, model.end, 0.0)
      # multi objective datawriter
        sampler.datawriter = DataWriter(prefix, model.params, model.begin, model.end, simthreshold_NS = 0.50, 
                                        simthreshold_pbias = 25.0, simthreshold_rsr = 0.70,
                                        with_valid_data = model.with_valid_data,
                                        shift_one_day = model.shift_one_day)
        # Now we can sample with the implemented Latin hypercube algorithm:
        sampler.sample(runs)
