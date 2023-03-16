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

def PID_RT(SP,PV,Man,MVMan,MVFF,Kc,Ti,Td,alpha,Ts,MVMin,MVMax,MV,MVP,MVI,MVD,E,ManFF=False,PVInit=0,method='EBD-EBD'):
    Tfd=alpha*Td
    if len(PV)==0:
        E.append(SP[-1]-PVInit)
    else:
        E.append (SP[-1]-PV[-1])
    #Vecteur MVP
    MVP.append(Kc*E[-1])

    #Vecteur MVI
    if len(MVI)==0:
        MVI.append((Kc*Ts/Ti)*E[-1]) #tjrs initialiser avec EBD
    else:
        if method=='TRAP-TRAP':
            MVI.append(MVI[-1]+(0.5*Kc*Ts/Ti)*(E[-1]+E[-2]))
        else:
            MVI.append(MVI[-1]+(Kc*Ts/Ti)*E[-1]) #EBD

    #Vecteur MVD
    if len(MVD)==0:
        MVD.append((Kc*Td/(Tfd+Ts))*(E[-1]-E[-2]))
    else:
        if method=='TRAP-TRAP':
            MVD.append((Tfd-(Ts/2))/(Tfd+(Ts/2))*MVD[-1]+(Kc*Td/(Tfd+(Ts/2)))*(E[-1]-E[-2]))
        else:
            MVD.append((Tfd/(Tfd+Ts))*MVD[-1]+(Kc*Td/(Tfd+Ts))*(E[-1]-E[-2]))

    #Vecteur MV
    MV.append(MVP+MVI+MVD)






 #--------------------------------------------------- ma version   
    if method=='TRAP-TRAP':
        MVP.append(Kc*E[-1])
        if len(MVI)==0:
            MVI.append(((Kc*Ts)/(2*Ti)))*(E[-1]+E[-2])
        else:
            MVI.append(MVI[-1]+((Kc*Ts)/(2*Ti)))*(E[-1]+E[-2])
        if len(MVD)==0:
            MVD.append(0)
        else:
           MVD.append( ((Tfd-(Ts/2))/(Tfd+Ts/2))*MVD[-1]+( (Kc*Td)/(Tfd+(Ts/2)) )*E[-1]-E[-2])

        
        if MVP[-1]+MVI[-1]+MVD[-1]>MVMax:
            MVI[-1]=MVMax-MVP[-1]-MVD[-1]
        elif MVP[-1]+MVI[-1]+MVD[-1]<MVMin:
            MVI[-1]=MVMin-MVP[-1]-MVD[-1]
        

        if Man[-1]:
            if ManFF:
                MVI[-1]=MVMan[-1]-MVP[-1]-MVD[-1]
            else:
                MVI[-1]=MVMan[-1]-MVP[-1]-MVD[-1]
    elif method=='EBD-EBD':
        MVP.append(Kc*E[-1])
        if len(MVI)==0:
            MVI.append(Kc*(Ts/Ti)*E[-1])
        else:
            MVI.append(MVI[-1]+Kc*(Ts/Ti)*E[-1])
        if len(MVD)==0:
            MVD.append(0)
        else:
            MVD.append((Tfd/(Tfd+Ts))*MVD[-1]+(Kc*Td/(Tfd+Ts))*(E[-1]-E[-2]))

        
        if MVP[-1]+MVI[-1]+MVD[-1]>MVMax:
            MVI[-1]=MVMax-MVP[-1]-MVD[-1]
        elif MVP[-1]+MVI[-1]+MVD[-1]<MVMin:
            MVI[-1]=MVMin-MVP[-1]-MVD[-1]
        

        if Man[-1]:
            if ManFF:
                MVI[-1]=MVMan[-1]-MVP[-1]-MVD[-1]
            else:
                MVI[-1]=MVMan[-1]-MVP[-1]-MVD[-1]



        MV.append(MVP[-1]+MVI[-1]+MVD[-1])