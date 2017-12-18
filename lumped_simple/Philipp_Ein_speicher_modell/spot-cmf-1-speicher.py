#!/usr/bin/env python
# coding: utf-8 

"""

"""

from __future__ import division
import datetime

import cmf
import spotpy
from spotpy.parameter import Uniform as param
import os
import numpy as np

# Kalibrierungs Zeitraum...
# 1979: Spin-Up
begin = 1980
end = 1985
# 1986 .. 1988: Validierung

# Anzahl der Modellläufe
# TODO: Anpassen an Rechenzeit
runs = 500 #1000


fnQ = 'GrebenauQTagMittel__1979_1990.txt'
fnT = 'Temp_max_min_avg_1979_1988.txt'
fnP = 'Prec_Grebenau_1979_1988.txt'


class SingleStorage(object):
    """
    Enthält das gesamte Modell
    """
    # TODO: Namen des Modells anpassen
    def __init__(self,begin,end):
        """
        Initialisiert das Modell und baut das Grundsetup zusammen
        
        MP111 - Änderungsbedarf: Sehr groß, hier wird es Euer Modell definiert
        Verständlichkeit: mittel

        """
        
        
        # TODO: Parameterliste erweitern und anpassen. 
        # Wichtig: Die Namen müssen genau die gleichen sein,
        # wie die Argumente in der Funktion setparameters.
        #
        # Parameter werden wie folgt definiert:
        # param(<Name>,<Minimum>,<Maximum>)
        self.params = [param('tr',0.,1000.),
                       param('ETV1',0.,1000.),
                       param('Vr',0.,1000.),
                       param('fETV0',0.,1.),
                       param('beta',0.3,5.0)
                       ]
        
        
        # Lädt die Treiber Daten    
        P,T,Tmin,Tmax,Q = self.loadPETQ()
        self.Q=Q
        # Nutze nur einen Kern - für unser Modell ist das schneller
        cmf.set_parallel_threads(1)
        
        # Erzeuge ein Projekt mit einer Zelle für lumped Modell
        self.project = cmf.project()
        p = self.project
        c = p.NewCell(0,0,0,1000)
        
        # Füge Speicher zur Zelle
        l=c.add_layer(1.0)
        # TODO: Weitere Layer / Speicher zur Zelle hinzufügen
        
        # Verdunstung
        cmf.HargreaveET(l,c.transpiration)

        # Outlet
        self.outlet = p.NewOutlet('outlet',10,0,0)
        
        # Regen
        self.makestations(P,T,Tmin,Tmax)
        self.project = p
        self.begin = begin
        self.end = end   
    
    def loadPETQ(self):
        """
        Lädt Klima und Abfluss aus den entsprechenden Dateien fnQ, fnT, fnP
        
        MP111 - Änderungsbedarf: Eigentlich keiner
        Verständlichkeit: gut

        Reference: http://fb09-pasig.umwelt.uni-giessen.de/cmf/wiki/CmfTutTestData
        
        """    
        # Fester Modell-Startpunkt
        begin = datetime.datetime(1979,1,1)
        step = datetime.timedelta(days=1)
        # Leere Zeitreihen 
        P = cmf.timeseries(begin, step)
        P.extend(float(Pstr) for Pstr in open(fnP))
        
        Q = cmf.timeseries(begin,step)
        Q.extend(float(Qstr) for Qstr in open(fnQ))
        # Convert m3/s to mm/day
        Q *= 86400 * 1e3 / (2976.41 * 1e6)
        T = cmf.timeseries(begin,step)
        Tmin = cmf.timeseries(begin,step)
        Tmax = cmf.timeseries(begin,step)
        
        # Durchlaufen der Zeilen in der Datei
        for line in open(fnT):
            columns = line.split('\t')
            if len(columns) == 3:
                Tmax.add(float(columns[0]))
                Tmin.add(float(columns[1]))
                T.add(float(columns[2]))
                
        return P,T,Tmin,Tmax,Q
        
    def makestations(self,P,T,Tmin,Tmax):
        """
        Erzeugt die Regenfall und Klimastation
        P = Zeitreihe Niederschlag
        T, Tmin, Tmax = Zeitreihe Tägliche Mitteltemperatur, Min und Max
        
        MP111 - Änderungsbedarf: Eigentlich keiner
        Verständlichkeit: gut
        
        Reference: http://fb09-pasig.umwelt.uni-giessen.de/cmf/wiki/CmfTutMeteostation        
        """
        rainstation = self.project.rainfall_stations.add('Grebenau avg',P,(0,0,0))
        self.project.use_nearest_rainfall()

        # Temperaturdaten
        meteo = self.project.meteo_stations.add_station('Grebenau avg',(0,0,0))
        meteo.T = T
        meteo.Tmin = Tmin
        meteo.Tmax = Tmax
        self.project.use_nearest_meteo()
        
        return rainstation

    def setparameters(self,tr,Vr=0.0,V0=1000.,beta=1.0,ETV1=500.,fETV0=0.0,initVol = 100.):
        """
        Setzt die Parameter, dabei werden parametrisierte Verbindungen neu erstellt
        
        MP111 - Änderungsbedarf: Sehr groß, hier werden alle Parameter des Modells gesetzt
        Verständlichkeit: mittel

        """
        # Ein paar Abkürzungen um Tipparbeit zu sparen
        c = self.project[0]
        outlet = self.outlet
        
        # Setze den Water-Uptakestress
        c.set_uptakestress(cmf.VolumeStress(ETV1,ETV1 * fETV0))
        cmf.kinematic_wave(c.layers[0],outlet,tr/V0,exponent=beta,residual=Vr,V0=V0)
        c.layers[0].volume = initVol
        # TODO: Alle weiteren Connections / Parameter von Eurem Modell aufbauen


    def runmodel(self,verbose=False):
        """
        Startet das Modell
        
        verbose = Wenn verbose = True, dann wird zu jedem Tag etwas ausgegeben
        
        MP111 - Änderungsbedarf: Gering, kann so bleiben, kann aber auch
                        um andere Ergebniszeitreihen ergänzt werden. Achtung,
                        falls ihr mehrere Outlets benutzt
        Verständlichkeit: gut
        
        Reference: http://fb09-pasig.umwelt.uni-giessen.de/cmf/wiki/CmfTutFluxes
        """
        # Erzeugt ein neues Lösungsverfahren
        solver = cmf.ImplicitEuler(self.project,1e-9)
        # Verkürzte schreibweise für die Zelle - spart tippen
        c = self.project[0]
        
        # Neue Zeitreihe für Modell-Ergebnisse (res - result)
        resQ = cmf.timeseries(self.begin,cmf.day)
        # Starte den solver und rechne in Tagesschritten
        for t in solver.run(self.project.meteo_stations[0].T.begin,self.end,cmf.day):
            # Fülle die Ergebnisse
            if t>=self.begin:
                resQ.add(self.outlet.waterbalance(t))
            # Gebe eine Statusmeldung auf den Bildschirm aus,
            # dann weiß man wo der solver gerade ist
            if verbose:
                print(t,'Q=%5.3f, P=%5.3f' % (resQ[t],c.get_rainfall(t)))
        # Gebe die gefüllten Ergebnis-Zeitreihen zurück
        return resQ
        
    def simulation(self,vector):
        """
        SpotPy erwartet eine Methode simulation. Diese methode ruft einfach
        setparameters und runmodel auf, so dass spotpy zufrieden ist
        
        MP111 - Änderungsbedarf: Keiner, außer bei runmodel ändert sich die return-Zeile
        Verständlichkeit: naja
        
        """        
        
        paramdict = dict((pp.name,v) for pp,v in zip(self.params,vector))
        self.setparameters(**paramdict)
        resQ = self.runmodel()
        return np.array(resQ)

    def evaluation(self):
        """
        Gehört auch zum spotpy-Interface. 

        MP111 - Änderungsbedarf: Keiner
        Verständlichkeit: Schlecht
        """
        return np.array(self.Q[self.begin:self.end + datetime.timedelta(days=1)])
    

    def parameters(self):
        """
        Gehört auch zum spotpy-Interface. 

        MP111 - Änderungsbedarf: Keiner
        Verständlichkeit: Schlecht
        """        
        return spotpy.parameter.generate(self.params)

    def objectivefunction(self,simulation,evaluation):
        """
        Gehört auch zum spotpy-Interface. 

        MP111 - Änderungsbedarf: Keiner
        Verständlichkeit: Mittel
        """        
        return spotpy.objectivefunctions.nashsutcliffe(evaluation,simulation)

    def plotsimulation(self,threshold):
        """
        Erzeugt einen plot der Behaviourial runs, des besten runs und
        der Beobachtung
        
        MP111 - Änderungsbedarf: Nicht sofort...
        Verständlichkeit: Mittel
        
        """
        #%%
        import pylab as plt
        simulations = np.loadtxt('lhs-1stor.csv',delimiter=',',skiprows=1)
        best_run = simulations[:,0].argmax()
        best = simulations[best_run,len(self.params)+1:-1]

        data = simulations[simulations[:,0]>=threshold, len(self.params)+1:-1]
       #%%
        p5 = np.percentile(data,5,0)
        p95 = np.percentile(data,95,0)
        x = plt.arange(len(best))
        plt.fill_between(x,p5,p95,facecolor='#ffff00',edgecolor='none')
        plt.plot(x,best,'r-')
        plt.plot(x,self.evaluation(),'k')
        plt.show()
#%%
        
        
        

# http://stackoverflow.com/questions/419163/what-does-if-name-main-do
if __name__ == '__main__':
    # Importiere Algorithmus
    from spotpy.algorithms import rope as Sampler

    # Finde heraus, ob das ganze parallel laufen soll (für Supercomputer)
    parallel = 'mpi' if 'OMPI_COMM_WORLD_SIZE' in os.environ else 'seq'
    
    # Create the spotted model    
    model = SingleStorage(datetime.datetime(begin,1,1),
                          datetime.datetime(end,12,31))
    if runs:
        sampler = Sampler(model, parallel=parallel, dbname='lhs-1stor', dbformat='csv', save_sim=True)
        # sampler.datawriter = DataWriter(model.params, model.begin, model.end, 0.0)
        # Now we can sample with the implemented Monte Carlo algortihm:
        sampler.sample(runs)
        
    # plottet das Ergebnis
    # TODO: Threshold anpassen für die graue Fläche im Plot. Eure Modelle schaffen hoffentlich deutlich mehr...
    model.plotsimulation(0.3)
    
