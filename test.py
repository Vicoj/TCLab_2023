
from signal import signal
from tokenize import Single
import numpy as np
from scipy.optimize import minimize
from matplotlib import colors as mcolors

import matplotlib.pyplot as plt

from matplotlib.widgets import Slider, Button, RadioButtons,TextBox,CheckButtons
import matplotlib.patches as mpatches
import datetime
import os
import tclab
import time

from xmlrpc.client import Boolean
import numpy as np

from IPython.display import display, clear_output

import package_Class
from package_Class import *





SIM = Simulation(1000,1,20,True,'CPL_PID_FF_MAN_MODE')
G = Graph(SIM)

SP = Path(SIM,{0:50, 250:50 , SIM.TSim: 50})
DV = Path(SIM,{0:30 ,500:70,SIM.TSim: 70})
MAN = Path(SIM,{0:1 ,50: 0,500: 0, SIM.TSim: 0})
MANV = Path(SIM,{0: 70,400:100,600:30, SIM.TSim: 100})

P = FirstOrder(SIM,0.6,200,1,50,SIM.PVInit)
D = FirstOrder(SIM,0.4,152,29,50,0)

FF = FeedForward(SIM,P,D,True)
PID = PID_Controller(SIM,5,200,0,0.14,[0,100],False,False)
#PID.IMC_tuning(P,0.1,'H')

if(SIM.sim == True):
    t = []
    for ti in SIM.t:
        t.append(ti)
        #Signals
        SP.RT(t)
        DV.RT(t)
        MAN.RT(t)
        MANV.RT(t)

        FF.RT(DV.Signal) # FeedForward
        PID.RT(SP.Signal,SIM.PV,MAN.Signal,MANV.Signal,FF.MVFF,'EBD-EBD')

        SIM.MV.append(PID.MVFB[-1]+FF.MVFF[-1])

        P.RT(SIM.MV,'EBD')
        D.RT(DV.Signal,'EBD')
        
        SIM.PV.append(P.PV[-1]+D.PV[-1])  # Point Value

if(SIM.sim == False):
    #Tc Lab
    LAB = tclab.TCLab()
    LABVal = LabValues(SIM,LAB)

    SIM.t = []
    start = time.time()
    delta = 0
    totalTime = 0
    last = time.time()

    while totalTime < SIM.TSim:
        if delta > SIM.Ts:
            last = time.time()
            SIM.t.append(round(totalTime,4))
            #Signals
            SP.RT(SIM.t)
            DV.RT(SIM.t)
            MAN.RT(SIM.t)
            MANV.RT(SIM.t)

            FF.RT(DV.Signal) # FeedForward
            PID.RT(SP.Signal,SIM.PV,MAN.Signal,MANV.Signal,[0],'EBD-EBD')
    
            SIM.MV.append(PID.MVFB[-1]+FF.MVFF[-1])

            LABVal.RT(SIM.MV,DV.Signal,D.point_fct)
            delta = 0
            SIM.updateBar()

        else :
            totalTime = time.time() - start
            delta = time.time() - last
    LAB.close()

# Signaux Afficher dans  Graph Binaires
SigValsBin = [
    Signal(MAN.Signal,'Manual Mode','-g'),
]

# Signaux Afficher dans  Graph 1 : Temperature
SigVals1 = [
    Signal(SP.Signal,'SP','-r'),
    Signal(SIM.PV,'PV','-b')
]

# Signaux Afficher dans Graph 2 : % de chauffe
SigVals2 = [
    Signal(SIM.MV,'MV','-g'),
    Signal(DV.Signal,'DV','-k'),
    Signal(PID.MVFB,'MVFB','-y'),
    Signal(FF.MVFF,'MVFF','-r'),
    Signal(PID.MVP,'MVP',':g'),
    Signal(PID.MVI,'MVI',':y'),
    Signal(PID.MVD,'MVD',':m'),
]

# Signaux enregistrer dans le .txt
SigSave = [
    Signal(SIM.MV,'MV','-b'),
    Signal(PID.MVP,'MVP','-b'),
    Signal(PID.MVI,'MVI','-b'),
    Signal(PID.MVD,'MVD','-b'),
    Signal(SP.Signal,'SP',':m'),
    Signal(SIM.PV,'PV',':m'),
    Signal(DV.Signal,'DV',':m'),

]

# Variables affichées sur le graph
varVals = [
    Variable(SIM.TSim,'Temps Sim [s]'),
    Variable(SIM.Ts,'Sampling [s]'),
    Variable(SIM.PVInit,'Pv Init [°C]'),

    Variable(PID.OLP,'Open Loop'),
    Variable(PID.ManFF,'Man FF'),

    Variable(PID.Kc,'Kc PID'),
    Variable(PID.Td,'Td PID'),
    Variable(PID.Ti,'Ti PID'),

    Variable(FF.active,'FF Enabled'),
    Variable(FF.T1p,'TLead P(s)'),
    Variable(FF.T2p,'TLag P(s)'),
    Variable(FF.T1d,'TLead D(s)'),
    Variable(FF.T2d,'TLag D(s)'),

    Variable(P.Kg,'P(s) Gain K'),
    #Variable(P.T,'P(s) Time T'),
    Variable(P.Theta,'P(s) Theta θ'),

    Variable(D.Kg,'D(s) Gain K'),
    #Variable(D.T,'D(s) Time T'),
    Variable(D.Theta,'D(s) Theta θ'),
]

G.show([SigVals1,SigVals2,SigSave],SigValsBin,varVals)
G.Bode(P,PID)