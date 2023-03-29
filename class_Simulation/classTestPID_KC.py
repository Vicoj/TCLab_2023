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
from package_Class import Simulation,Path,FirstOrder,SecondOrderPlusDelay,LeadLag,FeedForward,PID_Controller,Delay,Signal,Variable,Graph,LabValues


#Simulation Instance
SIM = Simulation(2000,1,26,True,'EXP_RESP_TO_DV_FF_AUTOM_2')
SIM2 = Simulation(2000,1,26,True,'EXP_RESP_TO_DV_FF_AUTOM_2')
SIM3 = Simulation(2000,1,26,True,'EXP_RESP_TO_DV_FF_AUTOM_2')# Graph Instance
G = Graph(SIM,'PID Control_')

# Path

SP = Path(SIM,{0: 50, 1000: 60, SIM.TSim: 60})
DV = Path(SIM,{0: 50, SIM.TSim: 50})
MAN = Path(SIM,{0: 0, SIM.TSim: 0})
MANV = Path(SIM,{0: 0, SIM.TSim: 0})
SP2 = Path(SIM2,{0: 50, 1000: 60, SIM.TSim: 60})
DV2 = Path(SIM2,{0: 50, SIM.TSim: 50})
MAN2= Path(SIM2,{0: 0, SIM.TSim: 0})
MANV2 = Path(SIM2,{0: 0, SIM.TSim: 0})

SP3 = Path(SIM3,{0: 50, 1000: 60, SIM.TSim: 60})
DV3 = Path(SIM3,{0: 50, SIM.TSim: 50})
MAN3 = Path(SIM3,{0: 0, SIM.TSim: 0})
MANV3 = Path(SIM3,{0: 0, SIM.TSim: 0})



# FO Process
P = FirstOrder(SIM,0.6522434279003099,245.9823790885576,0.649693920059717,50,SIM.PVInit)
P2 = FirstOrder(SIM2,0.6522434279003099,245.9823790885576,0.649693920059717,50,SIM.PVInit)
P3 = FirstOrder(SIM3,0.6522434279003099,245.9823790885576,0.649693920059717,50,SIM.PVInit)
D = FirstOrder(SIM,0.6156105636473335,387.0591022229922, 5.419428855220769,50,0)
D2 = FirstOrder(SIM2,0.6156105636473335,387.0591022229922, 5.419428855220769,50,0)
D3 = FirstOrder(SIM3,0.6156105636473335,387.0591022229922, 5.419428855220769,50,0)

# Feed Forward

FF = FeedForward(SIM,P,D,True)
PID = PID_Controller(SIM,3.8187,246.3072,0.8094,2,0,100,False,True)
PID2 = PID_Controller(SIM2,1,246.3072,0.8094,2,0,100,False,True)
PID3 = PID_Controller(SIM3,6,246.3072,0.8094,2,0,100,False,True)
#PID.IMC_tuning(P,0.4,'H')



if(SIM.sim == True):
    t = []
    for ti in SIM.t:
        t.append(ti)
        #Signals
        SP.RT(t)
        DV.RT(t)
        MAN.RT(t)
        MANV.RT(t)

        SP2.RT(t)
        DV2.RT(t)
        MAN2.RT(t)
        MANV2.RT(t)

        SP3.RT(t)
        DV3.RT(t)
        MAN3.RT(t)
        MANV3.RT(t)

        FF.RT(DV.Signal) # FeedForward
        
        PID.RT(SP.Signal,SIM.PV,MAN.Signal,MANV.Signal,FF.MVFF,'EBD-EBD')
        PID2.RT(SP2.Signal,SIM2.PV,MAN.Signal,MANV.Signal,FF.MVFF,'EBD-EBD')
        PID3.RT(SP3.Signal,SIM3.PV,MAN.Signal,MANV.Signal,FF.MVFF,'EBD-EBD')
        
        SIM.addMV(PID.MVFB,FF.MVFF) # Modified Value
        SIM2.addMV(PID2.MVFB,FF.MVFF)
        SIM3.addMV(PID3.MVFB,FF.MVFF)
        
        P.RT(SIM.MV,'EBD')
        D.RT(DV.Signal,'EBD')
        SIM.PV.append(P.PV[-1]+D.PV[-1]) # Point Value
        #SIM.updateBar()

        P2.RT(SIM2.MV,'EBD')
        D2.RT(DV.Signal,'EBD')
        SIM2.PV.append(P2.PV[-1]+D2.PV[-1]) # Point Value
        #SIM2.updateBar()

        P3.RT(SIM3.MV,'EBD')
        D3.RT(DV.Signal,'EBD')
        SIM3.PV.append(P3.PV[-1]+D3.PV[-1]) # Point Value
        #SIM3.updateBar()

    

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
            PID.RT(SP.Signal,SIM.PV,MAN.Signal,MANV.Signal,FF.MVFF,'EBD-EBD')

            SIM.addMV(PID.MVFB,FF.MVFF) # Modified Value

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
    Signal(SP.Signal,'SP','-m'),
    Signal(SIM.PV,'PV','-b'),
    Signal(SIM2.PV,'PV2','-r'),
    Signal(SIM3.PV,'PV3','-g'),
    #Signal(P.PV,'P(s)','--b'),
    #Signal(D.PV,'D(s)','--k'),
]
SigVals2 = [
    Signal(SIM.MV,'MV','-b'),
    Signal(SIM2.MV,'MV2','-r'),
    Signal(SIM3.MV,'MV3','-g'),
    #Signal(DV.Signal,'DV','-k'),
    #Signal(MANV.Signal,'MANVal','-m'),
    #Signal(FF.MVFF,'MVFF','-g'),
    #Signal(PID.MVFB,'MVFB','-y'),
    #Signal(PID.E,'E',':r'),
    #Signal(PID.MVP,'MVP',':g'),
    #Signal(PID.MVI,'MVI',':y'),
    #Signal(PID.MVD,'MVD',':m'),
    #Signal(DV.Signal,'DV','-k'),
]
SigSave = [
    Signal(SIM.MV,'MV','-b'),
    Signal(PID.MVP,'MVP',':b'),
    Signal(PID.MVI,'MVI',':y'),
    Signal(PID.MVD,'MVD',':m'),
    Signal(SP.Signal,'SP',':m'),
    Signal(SIM.PV,'PV',':m'),
    Signal(DV.Signal,'DV','-k'),
    Signal(MAN.Signal,'Man','-k'),
]
varVals = [
    Variable(SIM.TSim,'Temps Sim [s]'),
    Variable(SIM.Ts,'Sampling [s]'),
    Variable(SIM.PVInit,'Pv Init [°C]'),

    Variable(PID.OLP,'Open Loop'),
    Variable(PID.ManFF,'Man FF'),


    Variable(PID.Kc,'Kc PID'),
    Variable(PID.Td,'Td PID'),
    Variable(PID.Ti,'Ti PID'),
    Variable(PID.alpha,'alpha PID'),
    Variable(PID.gamma,'Gamma IMC'),

    Variable(FF.active,'FF Enabled'),
    Variable(FF.T1p,'TLead P(s)'),
    Variable(FF.T2p,'TLag P(s)'),
    Variable(FF.T1d,'TLead D(s)'),
    Variable(FF.T2d,'TLag D(s)'),

]

G.show([SigVals1,SigVals2,SigSave],SigValsBin,varVals)
#G.Bode(P,PID,'PID')
#G.Bode(P,PID,'P')
