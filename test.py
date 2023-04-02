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




SIM = Simulation(2000,1,20,False,'CPL_FF_MAN_MODE_ALL')
G = Graph(SIM)

SP = Path(SIM,{0:40, 500:90 , 1000:50 , SIM.TSim:50})
DV = Path(SIM,{0:70 , 750:30 , 1000:30 , 1500:50 , SIM.TSim: 50})
MAN = Path(SIM,{0:1 , 50:0 , 1750:1, SIM.TSim:1})
MANV = Path(SIM,{0: 70,400:100,600:30, SIM.TSim: 30})

'''
-MV
K: 0.6380496617368366
T1: 194.68042676186926
T2: 4.132530123359063e-11
theta: 5.142082869041576
-DV
K: 0.41729187227800746
T1: 165.26136274440154
T2: 17.712030277318842
theta: 7.447021274687229
'''

P = SecondOrder(SIM,0.64,194.7,4.1325e-11,5.14,50,SIM.PVInit)
D = SecondOrder(SIM,0.42,165.26,17.71,7.45,50,0)

FF = FeedForward(SIM,P,D,True)

PID = PID_Controller(SIM,5,100,1,0.14,[0,100],False,False)
PID.IMC_tuning(P,0.5,'B')



if(SIM.sim == True):
    t = []
    for ti in SIM.t:
        t.append(ti)
        #Signals
        SP.RT(t)
        DV.RT(t)
        MAN.RT(t)
        MANV.RT(t)

        FF.RT(DV.Signal)
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

            FF.RT(DV.Signal)
            PID.RT(SP.Signal,SIM.PV,MAN.Signal,MANV.Signal,FF.MVFF,'EBD-EBD')

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
    Signal(FF.MVFF,'MVFF','-m'),

    #Signal(PID.MVP,'MVP',':g'),
    #Signal(PID.MVI,'MVI',':y'),
    #Signal(PID.MVD,'MVD',':m'),
#
]

# Signaux enregistrer dans le .txt
SigSave = [
    Signal(SIM.MV,'MV','-b'),
    Signal(PID.MVP,'MVP','-b'),
    Signal(PID.MVI,'MVI','-b'),
    Signal(PID.MVD,'MVD','-b'),
    Signal(FF.MVFF,'MV_FF','-b'),
    Signal(SP.Signal,'SP',':m'),
    Signal(SIM.PV,'PV',':m'),
    Signal(DV.Signal,'DV',':m'),
    Signal(MAN.Signal,'Man','-k'),

]

# Variables affichées sur le graph
varVals = [
    Variable(SIM.TSim,'Temps Sim [s]'),
    Variable(SIM.Ts,'Sampling [s]'),
    Variable(SIM.PVInit,'Pv Init [°C]'),

    Variable(PID.OLP,'Open Loop'),
    Variable(PID.ManFF,'Man FF'),
    Variable(FF.active,'FF Enabled'),

    Variable(0,'IMC Case B'),
    Variable(PID.Kc,'Kc PID'),
    Variable(PID.Td,'Td PID'),
    Variable(PID.Ti,'Ti PID'),

    Variable(P.Kg,'P(s) Gain K'),
    Variable(P.T1,'P(s) Time T1'),
    Variable(P.T2,'P(s) Time T2'),
    Variable(P.Theta,'P(s) Theta θ'),

    Variable(D.Kg,'D(s) Gain K'),
    Variable(D.T1,'D(s) Time T1'),
    Variable(D.T2,'D(s) Time T2'),
    Variable(D.Theta,'D(s) Theta θ'),
]

G.show([SigVals1,SigVals2,SigSave],SigValsBin,varVals)