import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output

#-----------------------------------
def SelectPath_RT(path,time,signal):
    
    """
    The function "SelectPath_RT" needs to be included in a "for or while loop".
    
    :path: dictionary input describing a path in time. Example: path = {0: 0, 5: 1, 50: 2, 80: 3, 100: 3}
    :time: time vector.
    :signal: signal vector that is being constructed using the input "path" and the vector "time".
    
    The function "SelectPath_RT" takes the last element in the vector "time" and, given the input "path", it appends the correct value to the vector "signal".
    """    
    
    for timeKey in path:
        if(time[-1] >= timeKey):
            timeKeyPrevious = timeKey    
    
    value = path[timeKeyPrevious]
    signal.append(value)

#-----------------------------------
def Delay_RT(MV,theta,Ts,MV_Delay,MVInit=0):
    
    """
    The function "Delay_RT" needs to be included in a "for or while loop".
    
    :MV: input vector
    :theta: delay [s]
    :Ts: sampling period [s]
    :MV_Delay: delayed input vector
    :MVInit: (optional: default value is 0)
    
    The function "Delay_RT" appends a value to the vector "MV_Delay".
    The appended value corresponds to the value in the vector "MV" "theta" seconds ago.
    If "theta" is not a multiple of "Ts", "theta" is replaced by Ts*int(np.ceil(theta/Ts)), i.e. the closest multiple of "Ts" larger than "theta".
    If the value of the vector "input" "theta" seconds ago is not defined, the value "MVInit" is used.
    """
    
    NDelay = int(np.ceil(theta/Ts))
    if NDelay > len(MV)-1:
        MV_Delay.append(MVInit)
    else:    
        MV_Delay.append(MV[-NDelay-1])

#-----------------------------------     

def Leadlag(MV,Kp,Tlead,Tlag,Ts,PV,PVInit=0,method='EBD'):
    
    """
    The function "Leadlag" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain
    :Tlead: lead time constant [s]
    :Tlag:lag time constant [s]
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
    
    if (Tlag != 0):
        K = Ts/Tlag
        if len(PV) == 0:
            PV.append(PVInit)
        else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
            if method == 'EBD':
                PV.append((1/(1+K))*PV[-1] + (K*K/(1+K))*((1+Tlead/Ts)*MV[-1]- (Tlead/Ts)*MV[-2]))
            elif method == 'EFD':
                PV.append((1-K)*PV[-1] + (K*K)*((Tlead/Ts)*MV[-1] +(1-(Tlead/Ts))*MV[-2]))
            elif method == 'TRAP':
                pass          
            else:
                PV.append((1/(1+K))*PV[-1] + (K*K/(1+K))*((1+Tlead/Ts)*MV[-1]- (Tlead/Ts)*MV[-2]))
    else:
        PV.append(Kp*MV[-1])
#-----------------------------------     

def PID_RT(SP,PV,Man,MVMan,MVFF,Kc,Ti,Td,alpha,Ts,MVMin,MVMax,MV,MVP,MVI,MVD,E,ManFF=False,PVInit=0,method='TRAP-TRAP'):
    """"
    The function "PID_RT3 needs to be included in a "for or while loop"

    :SP: SP (or SetPOint) vector
    :PV: PV (or Process Value) vector
    :Man: Man (or Manual controller mode) vector (True or False)
    :MVMan: MVMan (or Manual value for MV) vector
    :MVFF: MVFF (or Feedforward) vector

    :KC: controller gain
    :Ti: integral time constant [s]
    :Td: Derivative time constant [s]
    :alpha: Tfd=alpha*Td where Thd is the derivative filter time constant[s]
    :Ts: sampling period [s]

    :MVMin: minimum value for MV (used for saturation and anti wind-up)
    :MVMax: maximum value for MV (used for saturation and anti wind-up)

    :MV: MV (or Manipulated Value) vector
    :MVP: MVP (or Proportional part of MV) vector
    :MVI: MVI (or Integral part of MV) vector
    :MVD: MVD (or Derivative part of MV) vector
    :E: E (or control Error) vector

    :ManFF: Activated FF in manual mode (optimal:default boolean value is False)
    :PVInit: Initial value for PV (optional:default value is 0): used if PID_RT is ran firts in the sequence and no value of PV is available yet.

    :method: discretisation method (optional:default value is 'EBD')
    EBD-EBD: EBD for integral action and EBD for derivative action
    EBD-TRAP: EBD for integral action and TRAP for derivative action
    TRAP-TRAP: TRAP for integral action and TRAP for derivative action
    TRAP-EBD: TRAP for integral action and EBD for derivative action

    The function "PID_RT" appends new values to the vectors "MV","MV","MVI" and"MVD".
    The appeded values are based on the PID algorithm, the controller mode, and feedforward.
    Note that saturationof "MV" within the limits[MVMin,MVMax] is implemented with andi wind-up
    """
    Tfd=alpha*Td
    if len(PV)==0:
        PV.append(PVInit)
        E.append(SP[-1]-PVInit)
    else:
        E.append (SP[-1]-PV[-1])
    #Vecteur MVP
    MVP.append(Kc*E[-1])

    #Vecteur MVI
    if len(MVI)==0:
        MVI.append((Kc*Ts/Ti)*E[-1]) #tjrs initialiser avec EBD
    else:
        if method.split("-")[0]=='TRAP':
            MVI.append(MVI[-1]+(0.5*Kc*Ts/Ti)*(E[-1]+E[-2]))
        else:
            MVI.append(MVI[-1]+(Kc*Ts/Ti)*E[-1]) #EBD

    #Vecteur MVD
    if len(MVD)==0:
        if len(E)==1:
            MVD.append((Kc*Td/(Tfd+Ts))*(E[-1]))
        else:
            MVD.append((Kc*Td/(Tfd+Ts))*(E[-1]-E[-2]))
    else:
        if len(E)==1:
            if method.split("-")[1]=='TRAP':
                MVD.append((Tfd-(Ts/2))/(Tfd+(Ts/2))*MVD[-1]+(Kc*Td/(Tfd+(Ts/2)))*(E[-1]))
            else:
                MVD.append((Tfd/(Tfd+Ts))*MVD[-1]+(Kc*Td/(Tfd+Ts))*(E[-1]))
        else:
            if method.split("-")[1]=='TRAP':
                MVD.append((Tfd-(Ts/2))/(Tfd+(Ts/2))*MVD[-1]+(Kc*Td/(Tfd+(Ts/2)))*(E[-1]-E[-2]))
            else:
                MVD.append((Tfd/(Tfd+Ts))*MVD[-1]+(Kc*Td/(Tfd+Ts))*(E[-1]-E[-2]))

    #calcul saturation, anti emballement, reset saturation integrateur

    #mode automatique
    if(not Man[-1]):
        #saturation
        if(MVP[-1]+MVI[-1]+MVD[-1] + MVFF[-1] < MVMin) :
            MVI[-1] = MVMin - MVP[-1] - MVD[-1] - MVFF[-1] #ecrasement valeur de MV
        
        elif (MVP[-1]+MVI[-1]+MVD[-1] + MVFF[-1]>=MVMax) :
            MVI[-1] = MVMax - MVP[-1] - MVD[-1] - MVFF[-1]
        MV.append(MVP[-1]+MVI[-1]+MVD[-1])

    #mode manuel
    else :
        if(not ManFF):
            MVI[-1]=MVMan[-1]-MVP[-1]-MVD[-1]-MVFF[-1]
        else:
            MVI[-1]=MVMan[-1]-MVP[-1]-MVD[-1]

        MV.append(MVP[-1]+MVI[-1]+MVD[-1])
    