from ast import Str
from cProfile import label
from turtle import color
from xmlrpc.client import Boolean
import numpy as np
import pandas as pd
from matplotlib import colors as mcolors
from matplotlib.widgets import Slider, Button, RadioButtons,TextBox,CheckButtons
import matplotlib.pyplot as plt
from IPython.display import display, clear_output
from datetime import datetime
import os
import tclab
import progressbar



class Simulation:
    def __init__(self,TSim,Ts,PVInit,sim:bool,name:str()):
        if (sim):
            self.name = 'SIM_' + name
        else:
            self.name = 'EXP_' + name
        self.TSim = TSim
        self.sim = sim
        self.Ts = Ts
        self.N = int(TSim/Ts) + 1 
        self.PVInit = PVInit
        self.PV = []
        self.MV = []

        self.t = self.calc_t()
        self.i = 0
        widgets = [' [',progressbar.Timer(format= 'elapsed time: %(elapsed)s'),'] ',
           progressbar.Bar('█'),' (',
           progressbar.ETA(), ') ',
          ]
        self.bar = progressbar.ProgressBar(max_value=self.N, 
                              widgets=widgets).start()

    def calc_t(self):
        t = []
        for i in range(0,self.N):
            t.append(i*self.Ts)

        return t

    def updateBar(self):
        self.bar.update(self.t[self.i])
        self.i += 1

    def addMV(self,MVFB,MVFF):
        
        self.MV.append(MVFB[-1]+MVFF[-1])
        return None

class Path:
    """
        The function "SelectPath_RT" needs to be included in a "for or while loop".

        :path: dictionary input describing a path in time. Example: path = {0: 0, 5: 1, 50: 2, 80: 3, 100: 3}
        :time: time vector.
        :signal: signal vector that is being constructed using the input "path" and the vector "time".

        The function "SelectPath_RT" takes the last element in the vector "time" and, given the input "path", it appends the correct value to the vector "signal".
    """    
    
    def __init__(self,S:Simulation,path):
        self.S = S
        self.path = path
        self.Signal = []

    def RT(self,time):

        """
        The function "SelectPath_RT" needs to be included in a "for or while loop".

        :path: dictionary input describing a path in time. Example: path = {0: 0, 5: 1, 50: 2, 80: 3, 100: 3}
        :time: time vector.
        :signal: signal vector that is being constructed using the input "path" and the vector "time".

        The function "SelectPath_RT" takes the last element in the vector "time" and, given the input "path", it appends the correct value to the vector "signal".
    """  

        for timeKey in self.path:
            if(time[-1] >= timeKey):
                timeKeyPrevious = timeKey    
    
        value = self.path[timeKeyPrevious]
        self.Signal.append(value)

class Delay:
    def __init__(self,S:Simulation,theta):
        self.S = S
        self.MVInit = 0 # Par Default
        self.theta = theta
        self.PV = []

    def RT(self,MV):

        """
        The function "Delay_RT" needs to be included in a "for or while loop".

        :MV: input vector
        :theta: delay [s]
        :Ts: sampling period [s]
        :MV_Delay: delayed input vector
        :MVInit: (optional: default value is 0)

        The function "Delay_RT" appends a value to the vector "MV_Delay".
        The appened value corresponds to the value in the vector "MV" "theta" seconds ago.
        If "theta" is not a multiple of "Ts", "theta" is replaced by Ts*int(np.ceil(theta/Ts)), i.e. the closest multiple of "Ts" larger than "theta".
        If the value of the vector "input" "theta" seconds ago is not defined, the value "MVInit" is used.
        """

        NDelay = int(np.ceil(self.theta/self.S.Ts))
        if NDelay > len(MV)-1:
            self.PV.append(self.MVInit)
        else:    
            self.PV.append(MV[-NDelay-1])

class FirstOrder:
    def __init__(self,S:Simulation,gain,Time,Theta_FO,point_fct,PVInit):
        self.S = S
        self.Kg = gain
        self.T = Time
        self.Theta = Theta_FO
        self.point_fct = point_fct
        self.PVInit = PVInit

        self.FO_delay = Delay(S,Theta_FO)



        self.PV = []

    def RT(self,MV,method):

        """
        The function "FO_RT" needs to be included in a "for or while loop".

        :MV: input vector
        :Kp: process gain
        :T: lag time constant [s]
        :Ts: sampling period [s]
        :PV: output vector
        :PVInit: (optional: default value is 0)
        :method: discretisation method (optional: default value is 'EBD')
            EBD: Euler Backward difference
            EFD: Euler Forward difference
            TRAP: Trapezoïdal method

        The function "FO_RT" appends a value to the output vector "PV".
        The appended value is obtained from a recurrent equation that depends on the discretisation method.
        """    
        self.FO_delay.RT(MV)

        if (self.T != 0):
            K = self.S.Ts/self.T
            if len(self.PV) == 0:
                self.PV.append(self.PVInit)
            else:
                if method == 'EBD':
                    Null = ((1/(1+K))*self.PV[-1] + (K*self.Kg/(1+K))*self.FO_delay.PV[-1])
                    self.PV.append((1/(1+K))*self.PV[-1] + (K*self.Kg/(1+K))*self.FO_delay.PV[-1])
                elif method == 'EFD':
                    self.PV.append((1-K)*self.PV[-1] + K*self.Kg*self.FO_delay.PV[-2])
                elif method == 'TRAP':
                    self.PV.append((1/(2*self.T+self.S.Ts))*((2*self.T-self.S.Ts)*self.PV[-1] + self.Kg*self.S.Ts*(self.FO_delay.PV[-1] + self.FO_delay.PV[-2])))            
                else:
                    self.PV.append((1/(1+K))*self.PV[-1] + (K*self.Kg/(1+K))*self.FO_delay.PV[-1])
        else:
            self.PV.append(self.Kg*self.FO_delay.PV[-1])

class SecondOrder:
    def __init__(self,S:Simulation,gain,T1,T2,Theta,point_fct,PVInit):

        self.Kg = gain
        self.T1 = T1
        self.T2 = T2
        self.Theta = Theta
        self.point_fct = point_fct
        self.PVInit = PVInit


        self.FO1 = FirstOrder(S,gain,T1,Theta,point_fct,PVInit)
        self.FO2 = FirstOrder(S,1,T2,0,point_fct,PVInit)
        self.SO_delay = Delay(S,Theta)

        self.PV = []

    def RT(self,MV,method):

        """
        The function "FO_RT" needs to be included in a "for or while loop".

        :MV: input vector
        :Kp: process gain
        :T: lag time constant [s]
        :Ts: sampling period [s]
        :PV: output vector
        :PVInit: (optional: default value is 0)
        :method: discretisation method (optional: default value is 'EBD')
            EBD: Euler Backward difference
            EFD: Euler Forward difference
            TRAP: Trapezoïdal method

        The function "FO_RT" appends a value to the output vector "PV".
        The appended value is obtained from a recurrent equation that depends on the discretisation method.
        """    
        #self.SO_delay.RT(MV)
        self.FO1.RT(MV,method)
        self.FO2.RT(self.FO1.PV,method)

        self.PV.append(self.FO2.PV[-1])
 
class LeadLag:
    def __init__(self,S:Simulation,K,TLead,TLag):
        self.S = S
        self.K = K
        self.TLead = TLead
        self.TLag = TLag

        self.PV = []

    def RT(self,MV,method):
        if (self.TLag != 0):
            K = self.S.Ts/self.TLag

            if len(self.PV) == 0:
                self.PV.append(self.S.PVInit)
            else: 
                if method == 'EBD':
                    self.PV.append((1/(1+K))*self.PV[-1] + (K*self.K/(1+K))*((1+self.TLead/self.S.Ts)*MV[-1]- (self.TLead/self.S.Ts)*MV[-2]))
                elif method == 'EFD':
                    self.PV.append((1-K)*self.PV[-1] + (K*self.K)*((self.TLead/self.S.Ts)*MV[-1] +(1-(self.TLead/self.S.Ts))*MV[-2]))
                elif method == 'TRAP':
                    pass
                else:
                    self.PV.append((1/(1+K))*self.PV[-1] + (K*self.K/(1+K))*((1+self.TLead/self.S.Ts)*MV[-1]- (self.TLead/self.S.Ts)*MV[-2]))
        else:
            self.PV.append(self.K*MV[-1])

class FeedForward:
    def __init__(self,S:Simulation,P,D,active:bool):
        self.active = active
        self.S = S
        self.P = P
        self.D = D

        self.Kd = D.Kg
        self.Kp = P.Kg

        if type(P) is FirstOrder :
            self.T1p = P.T
            self.T2p = 1
        else :
            self.T1p = P.T1
            self.T2p = P.T2

        if type(D) is FirstOrder :
            self.T1d = D.T
            self.T2d = 1
        else :
            self.T1d = P.T1
            self.T2d = P.T2
        

        self.DV0 = D.point_fct

        self.ThetaD = D.Theta
        self.ThetaP = P.Theta


        self.MVFF = []
        #Gain
        KFF1 = -(self.Kd/self.Kp) 
        KFF2 = active

        #Delay
        thetaFF = max(self.ThetaD-self.ThetaP,0)
        self.delayFF = Delay(S,thetaFF)
        
        #leadLag
        self.LL1 = LeadLag(S,KFF1,self.T1p,self.T1d)
        self.LL2 = LeadLag(S,KFF2,self.T2p,self.T2d)

    def RT(self,DV):
        
        #Dephasage
        PVFF = DV-self.DV0*np.ones_like(DV)
        self.delayFF.RT(PVFF)
    
        self.LL1.RT(self.delayFF.PV,'EBD')
        self.LL2.RT(self.LL1.PV,'EBD')
        
        
        if(self.active):
            self.MVFF.append(self.LL2.PV[-1])
        else :
            self.MVFF.append(0)

class PID_Controller:

    def __init__(self,S:Simulation,Kc,Ti,Td,alpha,MVMinMax:list,OLP:bool,ManFF:bool):
        
        self.S = S
        self.Kc = Kc
        self.Ti = Ti
        self.Td = Td
        self.alpha = alpha
        self.MVMin = MVMinMax[0]
        self.MVMax = MVMinMax[1]

        self.OLP = OLP
        self.ManFF = ManFF
        self.gamma = 0
        self.case = ''

        self.MVMan = []

        self.MVFB = []
        self.MVP = []
        self.MVI = []
        self.MVD = []
        self.E = []

    def IMC_tuning(self,P:FirstOrder, gamma, case:str()='G'):
    #theta process
    #Kp gain process
    #T1p = time constant process
    #gamma for desired closed loop time constant

        self.gamma = gamma
        self.case = case

        Tc = gamma*P.T # 0.2 <gamma< 0.9

        if case == "G" :
            self.Kc = P.T/(P.Kg*Tc+P.Theta)
            self.Ti = P.T
            self.Td = 0

        if case == "H" :
            self.Kc = (P.T+P.Theta/2)/(P.Kg*Tc+P.Theta/2)
            self.Ti = P.T+P.Theta/2
            self.Td = (P.T*P.Theta)/(2*Tc+P.Theta)


        

    def RT(self,SP,PV,MAN,MVMan,MVFF,method="EBD"):
        """
        The function "PID.RT" needs to be included in a "for or while loop".
        
        SP : Set Point
        PV : Processed Value
        MAN : Manual mode ON/OFF
        MVMan : Modified value set with Manual mode
        MVFF : Modified value calculated by the Feed Forward
        method : discretisation method (optional: default value is 'EBD')
        
        The function "PID.RT" appends a value to the following vectors :
        E, MVp, MVi, MVd, MVFB
        The appened values correspond to the values computed by a parallel PID controller.
        If MAN is OFF then the modified value is bounded by MVmin and MVmax. The MVi value is overwritten in order to respect the boundaries.
        If MAN is ON then the modified value is equal to MVMan. The MVi value is overwritten in order to avoid integrator wind-up.
        """
    #calcul de l'erreur SP-PV
    
        if(not self.OLP):
            if(len(PV)==0):
                self.E.append(SP[-1]-self.S.PVInit)
            else :
                self.E.append(SP[-1]-PV[-1])
        else:
            self.E.append(SP[-1])

        #calcul de MVp

        self.MVP.append(self.Kc*self.E[-1])

        #calcul MVi

        if(len(self.MVI)>0):
            self.MVI.append(self.MVI[-1]+(self.Kc*self.S.Ts*self.E[-1])/self.Ti)
        else :
            #self.MVI.append((self.Kc*(self.S.Ts/self.Ti)*(self.E[-1]*self.E[-2]))
            self.MVI.append(self.S.PVInit)

        #calcul MVd

        Tfd = self.alpha*self.Td
        if(self.Td>0):
            if(len(self.MVD)!=0):
                if(len(self.E)==1):
                    self.MVD.append((Tfd/(Tfd+self.S.Ts))*self.MVD[-1] + ((self.Kc*self.Td)/(Tfd+self.S.Ts))*(self.E[-1]))
                else:
                    self.MVD.append(( Tfd / (Tfd+self.S.Ts) )*self.MVD[-1] + ( (self.Kc*self.Td) / (Tfd+self.S.Ts) ) *(self.E[-1]-self.E[-2]))
            else :
                self.MVD.append(self.S.PVInit)

        if(self.Td ==0 ):
            self.MVD.append(0)

        #calcul saturation, anti emballement, reset saturation integrateur

        #mode automatique
        if(not MAN[-1]):
            #saturation
            if(self.MVP[-1]+self.MVI[-1]+self.MVD[-1] + MVFF[-1] < self.MVMin) :
                self.MVI[-1] = self.MVMin - self.MVP[-1] - self.MVD[-1] - MVFF[-1] #ecrasement valeur de MV
            
            elif (self.MVP[-1]+self.MVI[-1]+self.MVD[-1] + MVFF[-1]>=self.MVMax) :
                self.MVI[-1] = self.MVMax - self.MVP[-1] - self.MVD[-1] - MVFF[-1]
            self.MVFB.append(self.MVP[-1]+self.MVI[-1]+self.MVD[-1])

        #mode manuel
        else :
            if(not self.ManFF):
                self.MVI[-1]=MVMan[-1]-self.MVP[-1]-self.MVD[-1]-MVFF[-1]
            else:
                self.MVI[-1]=MVMan[-1]-self.MVP[-1]-self.MVD[-1]

            self.MVFB.append(self.MVP[-1]+self.MVI[-1]+self.MVD[-1])

class PD_Controller:
    def __init__(self,S:Simulation,Kc,Ti,Td,alpha,MVMin,MVMax,OLP,ManFF:bool):
        
        self.S = S
        self.Kc = Kc
        self.Ti = Ti
        self.Td = Td
        self.alpha = alpha
        self.MVMin = MVMin
        self.MVMax = MVMax
        self.OLP = OLP
        self.ManFF = ManFF
        self.gamma = 0
        self.case = ''

        self.MVMan = []

        self.MVFB = []
        self.MVP = []
        self.MVI = []
        self.MVD = []
        self.E = []

    def IMC_tuning(self,P:FirstOrder, gamma, case:str()):
    #theta process
    #Kp gain process
    #T1p = time constant process
    #gamma for desired closed loop time constant
        self.gamma = gamma
        self.case = case

        Tc = gamma*P.T # 0.2 <gamma< 0.9

        if case == "G" :
            self.Kc = P.T/(P.K*Tc+P.Theta)
            self.Ti = P.T+P.Theta/2
            self.Td = 0
        if case == "H" :
            self.Kc = (P.T+P.Theta/2)/(P.K*Tc+P.Theta/2)
            self.Ti = P.T+P.Theta/2
            self.Td = (P.T*P.Theta)/(2*Tc+P.Theta)


        

    def RT(self,SP,PV,MAN,MVMan,MVFF,method="EBD"):
        """
        The function "PID.RT" needs to be included in a "for or while loop".
        
        SP : Set Point
        PV : Processed Value
        MAN : Manual mode ON/OFF
        MVMan : Modified value set with Manual mode
        MVFF : Modified value calculated by the Feed Forward
        method : discretisation method (optional: default value is 'EBD')
        
        The function "PID.RT" appends a value to the following vectors :
        E, MVp, MVi, MVd, MVFB
        The appened values correspond to the values computed by a parallel PID controller.
        If MAN is OFF then the modified value is bounded by MVmin and MVmax. The MVi value is overwritten in order to respect the boundaries.
        If MAN is ON then the modified value is equal to MVMan. The MVi value is overwritten in order to avoid integrator wind-up.
        """
    #calcul de l'erreur SP-PV
    
        if(not self.OLP):
            if(len(PV)==0):
                self.E.append(SP[-1]-self.S.PVInit)
            else :
                self.E.append(SP[-1]-PV[-1])
        else:
            self.E.append(SP[-1])

        #calcul de MVp

        self.MVP.append(self.Kc*self.E[-1])

        #calcul MVd

        Tfd = self.alpha*self.Td
        if(self.Td>0):
            if(len(self.MVD)!=0):
                if(len(self.E)==1):
                    self.MVD.append((Tfd/(Tfd+self.S.Ts))*self.MVD[-1] + ((self.Kc*self.Td)/(Tfd+self.S.Ts))*(self.E[-1]))
                else:
                    self.MVD.append(( Tfd / (Tfd+self.S.Ts) )*self.MVD[-1] + ( (self.Kc*self.Td) / (Tfd+self.S.Ts) ) *(self.E[-1]-self.E[-2]))
            else :
                self.MVD.append(self.S.PVInit)

        #mode automatique
        if(not MAN[-1]):
            #saturation
            if(self.MVP[-1]+self.MVD[-1] + MVFF[-1] < self.MVMin) :
                self.MVFB.append(self.MVMin-MVFF[-1])
            
            elif (self.MVP[-1]+self.MVD[-1] + MVFF[-1]>=self.MVMax) :
                self.MVFB.append(self.MVMax -MVFF[-1])
            else :
                self.MVFB.append(self.MVP[-1]+self.MVD[-1])

        #mode manuel
        else :
            if(not self.ManFF):
                self.MVFB.append(MVMan[-1]-MVFF[-1])
            else:
                self.MVFB.append(MVMan[-1])

class P_Controller:

    def __init__(self,S:Simulation,Kc,MVMinMax:list,OLP:bool,ManFF:bool):
        
        self.S = S
        self.Kc = Kc

        self.MVMin = MVMinMax[0]
        self.MVMax = MVMinMax[1]

        self.OLP = OLP
        self.ManFF = ManFF

        self.MVMan = []

        self.MVFB = []
        self.MVP = []
        self.E = []


    def RT(self,SP,PV,method="EBD"):
        """
        The function "PID.RT" needs to be included in a "for or while loop".
        
        SP : Set Point
        PV : Processed Value
        MAN : Manual mode ON/OFF
        MVMan : Modified value set with Manual mode
        MVFF : Modified value calculated by the Feed Forward
        method : discretisation method (optional: default value is 'EBD')
        
        The function "PID.RT" appends a value to the following vectors :
        E, MVp, MVi, MVd, MVFB
        The appened values correspond to the values computed by a parallel PID controller.
        If MAN is OFF then the modified value is bounded by MVmin and MVmax. The MVi value is overwritten in order to respect the boundaries.
        If MAN is ON then the modified value is equal to MVMan. The MVi value is overwritten in order to avoid integrator wind-up.
        """
    #calcul de l'erreur SP-PV
    
        if(not self.OLP):
            if(len(PV)==0):
                self.E.append(SP[-1]-self.S.PVInit)
            else :
                self.E.append(SP[-1]-PV[-1])
        else:
            self.E.append(SP[-1])

        #calcul de MVp

        self.MVP.append(self.Kc*self.E[-1])

        if (self.MVP[-1] >  self.MVMax):
            self.MVP[-1] =  self.MVMax
        if (self.MVP[-1] <  self.MVMin):
            self.MVP[-1] =  self.MVMin
        self.MVFB.append(self.MVP[-1])

class LabValues:
    def __init__(self,S:Simulation,LAB:tclab):
        self.S = S
        self.LAB = LAB
        self.S.PVInit = self.LAB.T1
        self.Ps = []
        self.Ds = []


    def RT(self,MV,DV,DV0):

        self.LAB.Q1(MV[-1])
        self.LAB.Q2(DV[-1]-DV0)

        self.Ps.append(self.LAB.T1)
        self.Ds.append(self.LAB.T2-DV0)

        self.S.PV.append(self.Ps[-1])

class Signal:
    def __init__(self,Signal,name:Str(),color:Str()):
        self.Signal = Signal
        self.name = name
        self.color = color

class Variable:
    def __init__(self,var,name:Str()):
        self.var = var
        self.name = name

class Graph:
    def __init__(self,S:Simulation):

        self.S = S
        self.title = self.S.name
        
        

    def show(self,signals:list(),binSignals:list(),varVals:list()):
        self.fig, self.ax = plt.subplots(3, gridspec_kw={'height_ratios': [1,3,3]})
        self.signals = signals
        self.binSignals = binSignals
        self.varVals = varVals

        for bin in binSignals:
            self.ax[0].step(self.S.t,bin.Signal,bin.color,linewidth=2,label=bin.name,where='post')
        self.ax[0].set_ylabel('Valeurs Binaires (On/Off)')
        self.ax[0].set_xlabel('Temps [s]')
        self.ax[0].set_title(self.title)
        self.ax[0].legend(loc='best')


        for signal in signals[0]:
            self.ax1, =self.ax[1].step(self.S.t,signal.Signal,signal.color,linewidth=2,label=signal.name,where='post')
        self.ax[1].set_ylabel('Temperature [°C]')
        self.ax[1].set_xlabel('Temps [s]')
        #self.ax[1].set_ylim(-10,120)
        self.ax[1].legend(loc='best')

        for signal in signals[1]:
            self.ax2 =self.ax[2].step(self.S.t,signal.Signal,signal.color,linewidth=2,label=signal.name,where='post')

        self.ax[2].set_ylabel('Pourcentage de Chauffe [%]')
        self.ax[2].set_xlabel('Temps [s]')
        #self.ax[2].set_ylim([0, 100])
        self.ax[2].legend(loc='best')

        plt.subplots_adjust(left=0.05, bottom=0.05, right = 0.8,top=0.95,hspace=0.064)

        self.boxes = []
        i = 0.9
        for var in varVals:
            varbox = plt.axes([0.87,i , 0.05, 0.03])
            textVar =  TextBox(varbox, var.name+': ', initial=str(round(var.var,4)))
            i -= 0.04
            #textVar.on_submit(self.update)

            self.boxes.append(textVar)


    

        # Buttons
        saveax = plt.axes([0.93, 0.1, 0.05, 0.04]) #4-tuple of floats *rect* = [left, bottom, width, height]
        button_save = Button(saveax, 'Save', hovercolor='0.975')
        button_save.on_clicked(self.save)

        namebox = plt.axes([0.83, 0.1, 0.1, 0.04])
        self.text_box = TextBox(namebox, 'Name :', initial=self.S.name)

        closefig = plt.axes([0.83, 0.05, 0.1, 0.04])
        button_close = Button(closefig, 'Close', hovercolor='0.975')
        button_close.on_clicked(self.close)

        plt.get_current_fig_manager().full_screen_toggle()
        plt.show()        

    def saveFig(self,event):
        now = datetime.now()
        date_time = now.strftime("%Y-%m-%d-%Hh%M")
        self.fig.savefig("Output/Bode/Margin_Graph_"+self.text_boxMargin.text+ '_' + date_time)

    def save(self,event):

        
        now = datetime.now()
        date_time = now.strftime("%Y-%m-%d-%Hh%M")
        t = np.array(self.S.t)
        data = [t.T]
        data_names ='t,'

        for sig in self.signals[2]:
            arr = np.array(sig.Signal)
            data.append(arr.T)
            data_names += sig.name + ','
        
        my_data = np.vstack((data))
        my_data = my_data.T

        if (self.S.sim):
            self.fig.savefig("Output/SIM/PID_Graph_"+self.text_box.text+ '_' + date_time)
            nameFile = 'Data/SIM/PID_Graph_' + self.text_box.text + '_' + date_time + '.txt'
        else:
            self.fig.savefig("Output/EXP/PID_Graph_"+self.text_box.text+ '_' + date_time)
            nameFile = 'Data/EXP/PID_Graph_' + self.text_box.text + '_' + date_time + '.txt'

        if not os.path.exists('Data'):
            os.makedirs('Data')
        if not os.path.exists('Data/SIM'):
            os.makedirs('Data/SIM')
        if not os.path.exists('Data/EXP'):
            os.makedirs('Data/EXP')
        np.savetxt(nameFile,my_data,delimiter=',',header=data_names,comments='')
                   
    def close(self,event):
        plt.close()

    def Bode(self,P,PID:PID_Controller):
    
        """
        :P: Process as defined by the class "Process".
            Use the following command to define the default process which is simply a unit gain process:
                P = Process({})

            A delay, two lead time constants and 2 lag constants can be added.

            Use the following commands for a SOPDT process:
                P.parameters['Kp'] = 1.1
                P.parameters['Tlag1'] = 10.0
                P.parameters['Tlag2'] = 2.0
                P.parameters['theta'] = 2.0

            Use the following commands for a unit gain Lead-lag process:
                P.parameters['Tlag1'] = 10.0        
                P.parameters['Tlead1'] = 15.0        

        :omega: frequency vector (rad/s); generated by a command of the type "omega = np.logspace(-2, 2, 10000)".
        :Show: boolean value (optional: default value = True). If Show = True, the Bode diagram is shown. Otherwise Ps (P(j omega)) (vector of complex numbers) is returned.

        The function "Bode" generates the Bode diagram of the process P
        """  

        #Graduation Abssisse
        omega = np.logspace(-4, 1, 1000)
        s = 1j*omega

        Tfd = PID.alpha*PID.Td

        # PID TRSF
        Cs = PID.Kc * (1 + (1/(PID.Ti*s)) + (PID.Td*s)/(Tfd*s+1))
        
        if type(P) is FirstOrder:
            Ps = (P.Kg/(P.T*s + 1)) * np.exp(-P.Theta*s)
            courbe = 'First Order'
        if type(P) is SecondOrder:
            Ps = (P.Kg/((P.T1*s + 1) + (P.T2*s + 1))) * np.exp(-P.Theta*s)
            courbe = 'Second Order'

        


        Ls = np.multiply(Cs,Ps)

        Db0 = 0

        # On crée les valeurs ligne 0Db
        Gain_0Dby = np.array([Db0]*len(omega)) 

        
    
        self.fig, (ax_gain, ax_phase) = plt.subplots(2,1)

        # Gain part
        ax_gain.semilogx(omega,20*np.log10(np.abs(Ps)),label='P(s)')
        ax_gain.semilogx(omega,20*np.log10(np.abs(Ls)),label='L(s)')


        gain_min = np.min(20*np.log10(np.abs(Ls)/5))
        gain_max = np.max(20*np.log10(np.abs(Ls)*5))

        ax_gain.semilogx(omega,Gain_0Dby,label='0Db',color='k')
        ax_gain.semilogx(np.array([0]*len(Ps)),Gain_0Dby,label='Zero')
        
        ax_gain.set_xlim([np.min(omega), np.max(omega)])
        ax_gain.set_ylim([gain_min, gain_max])
        ax_gain.set_ylabel('Amplitude |P| [db]')
        ax_gain.set_title('Bode plot of Loop Gain')
        ax_gain.legend(loc='best')
        ax_gain.set_xscale("log")

        Phase_180y = np.array([-180]*len(Ps))

        # Phase part
        ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(Ps)),label='P(s)')
        ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(Ls)),label='Ls(s)')

        ax_phase.semilogx(omega,Phase_180y,label='-180°',color='k')

        ax_phase.set_xlim([np.min(omega), np.max(omega)])
        ph_min = np.min((180/np.pi)*np.unwrap(np.angle(Ls))) - 10
        ph_max = np.max((180/np.pi)*np.unwrap(np.angle(Ls))) + 10

        ax_phase.set_ylim([np.max([ph_min, -200]), ph_max])
        ax_phase.set_ylabel(r'Phase $\angle P$ [°]')
        ax_phase.legend(loc='best')


        #Find ωc
        val1 = 20*np.log10(np.abs(Ls))
        val2 = Gain_0Dby
        wcx = np.argwhere(np.diff(np.sign(val1 - val2))).flatten()

        #Find °c
        val3 = (180/np.pi)*np.unwrap(np.angle(Ls))
        val4 = Phase_180y
        tcx = np.argwhere(np.diff(np.sign(val3 - val4))).flatten()

        plt.subplots_adjust(left=0.05, bottom=0.05, right = 0.8,top=0.95,hspace=0.064)

        Kcbox = plt.axes([0.9,0.55 , 0.05, 0.03])
        textKc =  TextBox(Kcbox, 'Gain PID: ', initial=str(PID.Kc))

        Tibox = plt.axes([0.9,0.5 , 0.05, 0.03])
        textTi =  TextBox(Tibox, 'Ti PID: ', initial=str(PID.Ti))

        Tdbox = plt.axes([0.9,0.45 , 0.05, 0.03])
        textTd =  TextBox(Tdbox, 'Td PID: ', initial=str(PID.Td))

    
        if( len(wcx) == 0 and len(tcx) == 0 ):

            MGbox = plt.axes([0.9,0.65 , 0.05, 0.03])
            textMG =  TextBox(MGbox, 'Marge de gain [Db]: ', initial='None')

            AGbox = plt.axes([0.9,0.6 , 0.05, 0.03])
            textAG =  TextBox(AGbox, 'Marge de Phase [deg °]: ', initial='None')

            wcbox = plt.axes([0.9,0.75 , 0.05, 0.03])
            textWc =  TextBox(wcbox, 'Gain ωc: ', initial='None')

            acbox = plt.axes([0.9,0.7 , 0.05, 0.03])
            textAc =  TextBox(acbox, 'Angle ωc: ', initial='None')

        elif (len(wcx) == 0):

            AGbox = plt.axes([0.9,0.6 , 0.05, 0.03])
            textAG =  TextBox(AGbox, 'Marge de Phase [deg °]: ', initial='None')

            wcbox = plt.axes([0.9,0.75 , 0.05, 0.03])
            textWc =  TextBox(wcbox, 'Gain ωc: ', initial='None')

            ac = omega[tcx][0]
            GainMarg = val1[tcx][0]

            ax_gain.axvline(ac,color = 'y')
            ax_gain.scatter(ac,GainMarg,color = 'r')
            ax_gain.scatter(ac,Db0,color = 'r')

            ax_phase.axvline(ac,color = 'y')
            ax_phase.scatter(ac,-180,color = 'k')
            GainMarg = Db0 - GainMarg 

            acbox = plt.axes([0.9,0.7 , 0.05, 0.03])
            textAc =  TextBox(acbox, 'Angle ωc: ', initial=str(round(ac,4)))

            MGbox = plt.axes([0.9,0.65 , 0.05, 0.03])
            textMG =  TextBox(MGbox, 'Marge de gain [Db]: ', initial=str(round(GainMarg,4)))

        elif (len(tcx) == 0):

            acbox = plt.axes([0.9,0.7 , 0.05, 0.03])
            textAc =  TextBox(acbox, 'Angle ωc: ', initial='None')

            MGbox = plt.axes([0.9,0.65 , 0.05, 0.03])
            textMG =  TextBox(MGbox, 'Marge de gain [Db]: ', initial='None')

            wc = omega[wcx][0]
            phaseMarg = val3[wcx][0]

            ax_gain.axvline(wc,color = 'y')
            ax_gain.scatter(wc,Db0,color = 'k')

            ax_phase.axvline(wc,color = 'y')
            ax_phase.scatter(wc,-180,color = 'r')
            ax_phase.scatter(wc,phaseMarg,color = 'r')
            phaseMarg =  phaseMarg + 180

            wcbox = plt.axes([0.9,0.75 , 0.05, 0.03])
            textWc =  TextBox(wcbox, 'Gain ωc: ', initial=str(round(wc,4)))

            AGbox = plt.axes([0.9,0.6 , 0.05, 0.03])
            textAG =  TextBox(AGbox, 'Marge de Phase [deg °]: ', initial=str(round(phaseMarg,4)))
  
        else:
            ac = omega[tcx][0]
            GainMarg = val1[tcx][0]
            
            wc = omega[wcx][0]
            phaseMarg = val3[wcx][0]

            ax_gain.axvline(wc,color = 'y')
            ax_gain.axvline(ac,color = 'y')

            ax_gain.scatter(ac,GainMarg,color = 'r')
            ax_gain.scatter(ac,Db0,color = 'r')
            ax_gain.scatter(wc,Db0,color = 'k')

            ax_phase.axvline(wc,color = 'y')
            ax_phase.axvline(ac,color = 'y')

            ax_phase.scatter(ac,-180,color = 'k')
            ax_phase.scatter(wc,-180,color = 'r')
            ax_phase.scatter(wc,phaseMarg,color = 'r')


            GainMarg = Db0 - GainMarg 
            phaseMarg =  phaseMarg + 180
        

            wcbox = plt.axes([0.9,0.75 , 0.05, 0.03])
            textWc =  TextBox(wcbox, 'Gain ωc: ', initial=str(round(wc,4)))

            acbox = plt.axes([0.9,0.7 , 0.05, 0.03])
            textAc =  TextBox(acbox, 'Angle ωc: ', initial=str(round(ac,4)))

            MGbox = plt.axes([0.9,0.65 , 0.05, 0.03])
            textMG =  TextBox(MGbox, 'Marge de gain [Db]: ', initial=str(round(GainMarg,4)))

            AGbox = plt.axes([0.9,0.6 , 0.05, 0.03])
            textAG =  TextBox(AGbox, 'Marge de Phase [deg °]: ', initial=str(round(phaseMarg,4)))

        varbox = plt.axes([0.9,0.8 , 0.05, 0.03])
        textVar =  TextBox(varbox, 'Courbe: ', courbe)

        namebox = plt.axes([0.83, 0.15, 0.1, 0.04])
        self.text_boxMargin = TextBox(namebox,  'Name ', initial=self.title)

        saveax = plt.axes([0.93, 0.15, 0.05, 0.04]) #4-tuple of floats *rect* = [left, bottom, width, height]
        button_save = Button(saveax, 'Save', hovercolor='0.975')
        button_save.on_clicked(self.saveFig)

        closefig = plt.axes([0.83, 0.05, 0.1, 0.04])
        button_close = Button(closefig, 'Close', hovercolor='0.975')
        button_close.on_clicked(self.close)

        
        plt.get_current_fig_manager().full_screen_toggle()
        plt.show()

