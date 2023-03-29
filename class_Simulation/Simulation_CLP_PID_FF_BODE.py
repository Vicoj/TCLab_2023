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

SIM = Simulation(2000,1,20,True,'Exp_CLP_FF_IMC')

G = Graph(SIM,'PID Control_')

SP = Path(SIM,{0:40 , 1000:60 , SIM.TSim: 60})
DV = Path(SIM,{0:30 , 500: 30 , 1500:50 ,SIM.TSim: 50})
MAN = Path(SIM,{0:1 ,250: 0, SIM.TSim: 0})
MANV = Path(SIM,{0: 50, SIM.TSim: 30})

P = FirstOrder(SIM,0.63,146,1,50,SIM.PVInit)
D = FirstOrder(SIM,0.4,150,20,50,0)

P_delay = Delay(SIM,P.Theta)
D_delay = Delay(SIM,D.Theta)


FF = FeedForward(SIM,P,D,True)

PID = PID_Controller(SIM,50,1,5,0.14,[0,100],False,False)
#PID.IMC_tuning(P,0.5,'G')

G.Bode(P,PID,'Ls')