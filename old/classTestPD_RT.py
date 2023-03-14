from signal import signal
from tokenize import Single
import numpy as np
from scipy.optimize import minimize
from matplotlib import colors as mcolors
from matplotlib.widgets import Slider, Button, RadioButtons,TextBox,CheckButtons
import matplotlib.patches as mpatches
import datetime
import os
import tclab
import time


from xmlrpc.client import Boolean
import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output

import package_Class
from package_Class import Simulation,Path,FirstOrder,SecondOrderPlusDelay,LeadLag,FeedForward,PID_Controller,PD_Controller,Delay,Signal,Variable,Graph,LabValues


#Simulation Instance
SIM = Simulation(2000,1,26,False,'EXP_STEP_AUTO')

# Graph Instance
G = Graph(SIM,'PD_Control_')

# Path
SP = Path(SIM,{0: 50, 1000 : 60, SIM.TSim: 50})
DV = Path(SIM,{0: 50, SIM.TSim: 60})
MAN = Path(SIM,{0: 0, SIM.TSim: 0})
MANV = Path(SIM,{0: 50, SIM.TSim: 50})


# FO Process
P = FirstOrder(SIM,0.6522434279003099,245.9823790885576,0.649693920059717,50,SIM.PVInit)
D = FirstOrder(SIM,0.6156105636473335,387.0591022229922, 5.419428855220769,50,0)

# Feed Forward

FF = FeedForward(SIM,P,D,True)
PD = PD_Controller(SIM,1.69,141,5,2,0,100,False,True)

PD.IMC_tuning(P,0.4,'H')


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
        PD.RT(SP.Signal,SIM.PV,MAN.Signal,MANV.Signal,FF.MVFF,'EBD-EBD')

        SIM.addMV(PD.MVFB,FF.MVFF) # Modified Value

        P.RT(SIM.MV,'EBD')
        D.RT(DV.Signal,'EBD')
        SIM.PV.append(P.PV[-1]+D.PV[-1]) # Point Value
        SIM.updateBar()

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
            PD.RT(SP.Signal,SIM.PV,MAN.Signal,MANV.Signal,FF.MVFF,'EBD-EBD')

            SIM.addMV(PD.MVFB,FF.MVFF) # Modified Value

            LABVal.RT(SIM.MV,DV.Signal,D.point_fct)
            delta = 0
            SIM.updateBar()

        else :
            totalTime = time.time() - start
            delta = time.time() - last
    LAB.close()
        
SigValsBin = [
    Signal(MAN.Signal,'Manual Mode','-g')
]
SigVals1 = [
    Signal(SP.Signal,'SP','-r'),
    Signal(SIM.PV,'PV','-b'),
    #Signal(P.PV,'P(s)','--b'),
    #Signal(D.PV,'D(s)','--k'),
]
SigVals2 = [
    Signal(SIM.MV,'MV','-b'),
    Signal(DV.Signal,'DV','-k'),
    #Signal(MANV.Signal,'MANVal','-m'),
    #Signal(FF.MVFF,'MVFF','-g'),
    #Signal(PID.MVFB,'MVFB','-y'),
    #Signal(PID.E,'E',':r'),
    Signal(PD.MVP,'MVP',':b'),
    #Signal(PID.MVI,'MVI',':y'),
    Signal(PD.MVD,'MVD',':m'),
    #Signal(DV.Signal,'DV','-k'),
]
SigSave = [
    Signal(SIM.MV,'MV','-b'),
    Signal(PD.MVP,'MVP',':b'),
    Signal(PD.MVD,'MVD',':m'),
    Signal(SP.Signal,'SP',':m'),
    Signal(SIM.PV,'PV',':m'),
    Signal(DV.Signal,'DV','-k'),
    Signal(MAN.Signal,'Man','-k'),
]
varVals = [
    Variable(SIM.TSim,'Temps Sim [s]'),
    Variable(SIM.Ts,'Sampling [s]'),
    Variable(SIM.PVInit,'Pv Init [°C]'),

    Variable(PD.OLP,'Open Loop'),
    Variable(PD.ManFF,'Man FF'),


    Variable(PD.Kc,'Kc PID'),
    Variable(PD.Td,'Td PID'),
    Variable(PD.Ti,'Ti PID'),
    Variable(PD.gamma,'Gamma IMC'),

    Variable(FF.active,'FF Enabled'),
    Variable(FF.T1p,'TLead P(s)'),
    Variable(FF.T2p,'TLag P(s)'),
    Variable(FF.T1d,'TLead D(s)'),
    Variable(FF.T2d,'TLag D(s)'),

]

G.show([SigVals1,SigVals2,SigSave],SigValsBin,varVals)
#G.Bode(P,PID,'PID')
#G.Bode(P,PID,'P')
