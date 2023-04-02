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
        else: 
            if method == 'EBD':
                PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*((1+Tlead/Ts)*MV[-1]- (Tlead/Ts)*MV[-2]))
            elif method == 'EFD':
                PV.append((1-K)*PV[-1] + (K*Kp)*((Tlead/Ts)*MV[-1] +(1-(Tlead/Ts))*MV[-2]))
            elif method == 'TRAP':
                pass  
            else:
                PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*((1+Tlead/Ts)*MV[-1]- (Tlead/Ts)*MV[-2]))
    else:
        PV.append(Kp*MV[-1])
#-----------------------------------     

def PID_RT(SP,PV,Man,MVMan,MVFF,Kc,Ti,Td,alpha,Ts,MVMin,MVMax,MV,MVP,MVI,MVD,E,ManFF=False,PVInit=0,method='TRAP-TRAP'):
    """"
    The function "PID_RT needs to be included in a "for or while loop"

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




    #mode automatique
    if(not Man[-1]):
        #saturation
        if(MVP[-1]+MVI[-1]+MVD[-1] + MVFF[-1] < MVMin) :
            MVI[-1] = MVMin - MVP[-1] - MVD[-1] - MVFF[-1] #ecrasement valeur de MV
        
        elif (MVP[-1]+MVI[-1]+MVD[-1] + MVFF[-1]>=MVMax) :
            MVI[-1] = MVMax - MVP[-1] - MVD[-1] - MVFF[-1]
        MV.append(MVP[-1]+MVI[-1]+MVD[-1]+MVFF[-1])

    #mode manuel
    else :
        if(not ManFF):
            MVI[-1]=MVMan[-1]-MVP[-1]-MVD[-1]-MVFF[-1]
        else:
            MVI[-1]=MVMan[-1]-MVP[-1]-MVD[-1]

        MV.append(MVP[-1]+MVI[-1]+MVD[-1])


#----------------------------------- 

def IMC_TUNNING(Kp,theta,T1,T2,gamma,method='SOPDT'):

    """
    parameters :
    : Kp : process gain
    : theta : process offset
    : T1 : processus first time constant
    : T2 : processus second time constant
    : gamma : constant used to get the closed loop time constant
    : method : model of the IMC_Tunning (optional:default value is 'SOPDT')
    returns Kc, Ti, Td
    """

    Tc= gamma*T1 #Tclp
    if method=='FOPDT':
        Ti=T1
        Td=0
        Kc=(Ti/(Tc+theta)*Kp)
        return (Kc,Ti,Td)
    elif method=='SOPDT':
        Ti=T1+(theta/2)
        Td=(T1*theta)/(2*T1+theta)
        Kc=(T1+(theta/2))/(Tc+(theta/2)*Kp)
        return(Kc,Ti,Td)

#----------------------------------- 
     
class Controller:
    
    def __init__(self, parameters):
        
        self.parameters = parameters
        self.parameters['Kc'] = parameters['Kc'] if 'Kc' in parameters else 1.0
        self.parameters['alpha'] = parameters['alpha'] if 'alpha' in parameters else 0.0
        self.parameters['Ti'] = parameters['Ti'] if 'Ti' in parameters else 0.0
        self.parameters['Td'] = parameters['Td'] if 'Td' in parameters else 0.0

#-----------------------------------   

def Margins(P,C,omega, Show = True):
    
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
    
    :C: Process as defined by the class "Controller".
        Use the following command to define the default controller:
            C = Controller({})
        Use the following commands for a PID process:
            C.parameters['Kc'] = 1.1
            C.parameters['Ti'] = 10.0
            C.parameters['Td'] = 2.0
            C.parameters['alpha'] = 2.0      
        
    :omega: frequency vector (rad/s); generated by a command of the type "omega = np.logspace(-2, 2, 10000)".
    :Show: boolean value (optional: default value = True). If Show = True, the Bode diagram is shown. Otherwise Ps (P(j omega)) (vector of complex numbers) is returned.
    
    The function "Bode" generates the Bode diagram of the process P
    """     
    
    s = 1j*omega
    
    Ptheta = np.exp(-P.parameters['theta']*s)
    PGain = P.parameters['Kp']*np.ones_like(Ptheta)
    PLag1 = 1/(P.parameters['Tlag1']*s + 1)
    PLag2 = 1/(P.parameters['Tlag2']*s + 1)
    PLead1 = P.parameters['Tlead1']*s + 1
    PLead2 = P.parameters['Tlead2']*s + 1
    
    Ps = np.multiply(Ptheta,PGain)
    Ps = np.multiply(Ps,PLag1)
    Ps = np.multiply(Ps,PLag2)
    Ps = np.multiply(Ps,PLead1)
    Ps = np.multiply(Ps,PLead2)

    Cs = C.parameters['Kc']*(1 + (1/C.parameters['Ti']*s) + (C.parameters['Td']*s)/(C.parameters['alpha'] * C.parameters['Td'] * s + 1))
    
    Ls=np.multiply(Ps,Cs)

    if Show == True:
    
        fig, (ax_gain, ax_phase) = plt.subplots(2,1)
        fig.set_figheight(12)
        fig.set_figwidth(22)

        gain=20*np.log10(np.abs(Ls))
        phase=(180/np.pi)*np.unwrap(np.angle(Ls))
        gain_0 = np.argmin(np.abs(gain))
        phase_180 = np.argmin(np.abs(phase + 180))
        gain_interserction=omega[gain_0]
        phase_intersectoin=omega[phase_180]

        
        # Gain part
        pt1=[phase_intersectoin,gain[phase_180]]
        pt2=[phase_intersectoin,0]
        xgain_value = [pt1[0], pt2[0]]
        ygain_value = [pt1[1], pt2[1]]
        
        ax_gain.plot(xgain_value, ygain_value, color='orange')
        ax_gain.axhline(y=0, color='red',label='0dB')
        ax_gain.semilogx(omega,gain,label='L(s)')
        """ax_gain.semilogx(omega,20*np.log10(np.abs(PGain)),label='Pgain')
        if P.parameters['theta'] > 0:
            ax_gain.semilogx(omega,20*np.log10(np.abs(Ptheta)),label='Ptheta(s)')
        if P.parameters['Tlag1'] > 0:
            ax_gain.semilogx(omega,20*np.log10(np.abs(PLag1)),label='PLag1(s)')
        if P.parameters['Tlag2'] > 0:        
            ax_gain.semilogx(omega,20*np.log10(np.abs(PLag2)),label='PLag2(s)')
        if P.parameters['Tlead1'] > 0:        
            ax_gain.semilogx(omega,20*np.log10(np.abs(PLead1)),label='PLead1(s)')
        if P.parameters['Tlead2'] > 0:    
            ax_gain.semilogx(omega,20*np.log10(np.abs(PLead2)),label='PLead2(s)')    """
        gain_min = np.min(20*np.log10(np.abs(Ls)/5))
        gain_max = np.max(20*np.log10(np.abs(Ls)*5))
        ax_gain.set_xlim([np.min(omega), np.max(omega)])
        ax_gain.set_ylim([gain_min, gain_max])
        ax_gain.set_ylabel('Amplitude |L| [db]')
        ax_gain.set_title('Bode plot of L')
        ax_gain.legend(loc='best')
    
        # Phase part
        pt1 = [gain_interserction, -180]
        pt2 = [gain_interserction, phase[gain_0]]
        xphase_value = [pt1[0], pt2[0]]
        yphase_value = [pt1[1], pt2[1]]

        ax_phase.plot(xphase_value, yphase_value, color='orange')
        ax_phase.semilogx(omega, phase,label='L(s)')
        ax_phase.axhline(y=-180, color='red', label='-180°')

        """ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(PGain)),label='Pgain')
        if P.parameters['theta'] > 0:    
            ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(Ptheta)),label='Ptheta(s)')
        if P.parameters['Tlag1'] > 0:        
            ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(PLag1)),label='PLag1(s)')
        if P.parameters['Tlag2'] > 0:        
            ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(PLag2)),label='PLag2(s)')
        if P.parameters['Tlead1'] > 0:        
            ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(PLead1)),label='PLead1(s)')
        if P.parameters['Tlead2'] > 0:        
            ax_phase.semilogx(omega, (180/np.pi)*np.unwrap(np.angle(PLead2)),label='PLead2(s)')"""    
        ax_phase.set_xlim([np.min(omega), np.max(omega)])
        ph_min = np.min((180/np.pi)*np.unwrap(np.angle(Ls))) - 10
        ph_max = np.max((180/np.pi)*np.unwrap(np.angle(Ls))) + 10
        ax_phase.set_ylim([np.max([ph_min, -200]), ph_max])
        ax_phase.set_ylabel(r'Phase $\angle P$ [°]')
        ax_phase.legend(loc='best')
        print("gain: ",np.around(-gain[phase_180])," dB")
        print("phase: ",np.around(180+phase[gain_0]),"°")
   
    else:
        
        
        return Ls
    