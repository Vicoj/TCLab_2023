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
from package_Class import Simulation,Path,FirstOrder,SecondOrderPlusDelay,LeadLag,FeedForward,PID_Controller,Delay,Signal,Variable,Graph,LabValues

SIM = Simulation(2000,1,20,True,'Experience_OPL_FF_NO_PID')
G = Graph(SIM,'')

SP = Path(SIM,{0:30 , 1000:70 , SIM.TSim: 70})
DV = Path(SIM,{0: 50 ,1500:70 ,SIM.TSim: 70})

P = FirstOrder(SIM,0.6576982674557728,201.32536449625266,10,50,SIM.PVInit)
D = FirstOrder(SIM,0.4022800628236135,151.76786567122514,100,50,0)

P_delay = Delay(SIM,P.Theta)
D_delay = Delay(SIM,D.Theta)

FF = FeedForward(SIM,P,D,True)


if(SIM.sim == True):
    t = []
    for ti in SIM.t:
        t.append(ti)
        #Signals
        SP.RT(t)    # On crée le Set point
        DV.RT(t)    # On crée le Disturbance Load

        FF.RT(DV.Signal) # FeedForward

        SIM.addMV(SP.Signal,FF.MVFF) # Modified Value

        P_delay.RT(SIM.MV)
        D_delay.RT(DV.Signal)

        P.RT(P_delay.PV,'EBD') # P.RT(SIM.MV,'EBD')
        D.RT(D_delay.PV,'EBD')
 
        SIM.PV.append(P.PV[-1]+D.PV[-1]) # Point Value
        SIM.updateBar()

# Signaux Afficher dans  Graph Binaires
SigValsBin = [
]

# Signaux Afficher dans  Graph 1 : Temperature
SigVals1 = [
    Signal(SP.Signal,'SP','-r'),
    Signal(SIM.PV,'PV','-b'),
    Signal(DV.Signal,'DV','-k'),

    Signal(P.PV,'P(s)','--b'),
    Signal(D.PV,'D(s)','--k'),
]

# Signaux Afficher dans Graph 2 : % de chauffe
SigVals2 = [
    Signal(SIM.MV,'MV','-b'),
    Signal(DV.Signal,'DV','-k'),
    #Signal(MANV.Signal,'MANVal','-m'),

    Signal(FF.MVFF,'MVFF','-g'),
    Signal(FF.LL1.PV,'FF_LL1','--b'),
    Signal(FF.LL2.PV,'FF_LL2','--r'),
    Signal(FF.delayFF.PV,'FF_delay','--m'),

    #Signal(PID.MVFB,'MVFB','-y'),
    #Signal(PID.E,'E',':r'),
    #Signal(PID.MVP,'MVP',':g'),
    #Signal(PID.MVI,'MVI',':y'),
    #Signal(PID.MVD,'MVD',':m'),
]

# Signaux enregistrer dans le .txt
SigSave = [
    Signal(SIM.MV,'MV','-b'),
    #Signal(PID.MVP,'MVP',':b'),
    #Signal(PID.MVI,'MVI',':y'),
    #Signal(PID.MVD,'MVD',':m'),
    Signal(SP.Signal,'SP',':m'),
    Signal(SIM.PV,'PV',':m'),
    Signal(DV.Signal,'DV','-k'),
]

# Variables affichées sur le graph
varVals = [
    Variable(SIM.TSim,'Temps Sim [s]'),
    Variable(SIM.Ts,'Sampling [s]'),
    Variable(SIM.PVInit,'Pv Init [°C]'),

    #Variable(PID.OLP,'Open Loop'),
    #Variable(PID.ManFF,'Man FF'),
#
#
    #Variable(PID.Kc,'Kc PID'),
    #Variable(PID.Td,'Td PID'),
    #Variable(PID.Ti,'Ti PID'),
    #Variable(PID.gamma,'Gamma IMC'),
#
    Variable(FF.active,'FF Enabled'),
    Variable(FF.T1p,'TLead P(s)'),
    Variable(FF.T2p,'TLag P(s)'),
    Variable(FF.T1d,'TLead D(s)'),
    Variable(FF.T2d,'TLag D(s)'),
    Variable(FF.delayFF.theta,'Theta_FFD'),

]


G.show([SigVals1,SigVals2,SigSave],SigValsBin,varVals)
