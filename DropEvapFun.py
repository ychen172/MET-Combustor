import numpy as np
import cantera as ct
from scipy.optimize import fsolve
from scipy.integrate import odeint
#Compute Forced Convection Term
def DelTMCalc(MachLiner,Gamma,Tinf,Pref,YOxi,Ru,rd,gas):
    Tstatic = Tinf/(1 + 0.5*(Gamma-1)*(MachLiner**2))
    gas.TPY = Tinf,Pref,YOxi
    MWOxi = gas.mean_molecular_weight/1000 #kg/mole
    ROxi = Ru/MWOxi #Gas constant of fuel J/kg/K
    velRel = np.sqrt(Gamma*ROxi*Tstatic)*MachLiner #Approximated Relative Velocity
    Re = (gas.density_mass*velRel*(2*rd))/gas.viscosity #Relative Reynold Number
    Pr = (gas.viscosity*gas.cp_mass)/gas.thermal_conductivity #Relative Prandtl Number
    NuSh = 2+(0.555*Re**(1/2)*Pr**(1/3))/((1+1.232/(Re*Pr**(4/3)))**(1/2))
    DelTM = (NuSh/(NuSh-2))*rd #Unity Lewis Number Nu = Sh and so DelT = DelM
    return DelTM, NuSh, Re, Pr, velRel

def Objective(vars,extArgs):
    YOxi    = extArgs[0] #Composition of Oxidizer
    YFue    = extArgs[1] #Composition of Fuel
    gas     = extArgs[2]
    Pref    = extArgs[3] #Constant Pressure
    Tinf    = extArgs[4] #Temperature Farfield
    Tdrop   = extArgs[5] #Droplet Inner Temperature
    YOxiInf = extArgs[6] #Mass Fraction of Oxidizer in FarField
    rd      = extArgs[7] #Droplet Radius
    DelTM   = extArgs[8] #FarField Radius
    cpl     = extArgs[9] #Liquid specific heat
    LHVapor = extArgs[10] #Latent heat of vaporization
    LHVheat = extArgs[11] #Lower heating value
    FAst    = extArgs[12] #Stoichiometric Fuel Air Ratio
    TsatRef = extArgs[13] #Reference saturation temperature
    PsatRef = extArgs[14] #Reference saturation pressure
    MWFue   = extArgs[15] #Fuel Molecular Weight
    RFue    = extArgs[16] #Fuel Gas Constant
    MWPro   = extArgs[17] #Product Molecular Weight

    mdotF = vars[0] #Mass Flow Rate of Fuel
    rf    = vars[1] #Radius of Flame Sheet
    Tf    = vars[2] #Flamelet Temperature
    Ts    = vars[3] #Droplet Surface Temperature
    YFs   = vars[4] #Droplet Surface MassFraction

    #Average Property
    Tave = np.max([0.5*(Ts+Tf),Tdrop]) #Law and Williams (Prevent extreme low temperature)
    gas.TPY = Tave,Pref,YFue
    cpg  = gas.cp_mass #Law and Williams
    kF = gas.thermal_conductivity
    gas.TPY = Tave,Pref,YOxi
    kOx = gas.thermal_conductivity 
    kg = 0.4*kF + 0.6*kOx #Law and Williams
    rhoDg = kg/cpg # Unity Lewis Number
    #Solve Five Equ Five Unknown
    ZF = 1/(4*np.pi*rhoDg)
    ZT = ZF #Unity Lewis Number
    #Mass Residuals
    ResInnerMass = 1-(np.exp(-ZF*mdotF/rd)/np.exp(-ZF*mdotF/rf))  -   YFs
    ResOuterMass = (1/FAst)*(np.exp(-mdotF*ZF/DelTM)/np.exp(-mdotF*ZF/rf) - 1)  -  YOxiInf
    #Inner Energy Equation
    qil = cpl*(Ts-Tdrop) 
    ResInnerEnergy = (cpg*(Tf-Ts)*np.exp(-ZT*mdotF/rd))/((qil+LHVapor)*(np.exp(-ZT*mdotF/rd)-np.exp(-ZT*mdotF/rf)))  +  1
    #Outer Energy Equation
    TGradIn = ((Ts-Tf)*np.exp(-ZT*mdotF/rf))/(np.exp(-ZT*mdotF/rd)-np.exp(-ZT*mdotF/rf))
    TGradOut = ((Tf-Tinf)*np.exp(-ZT*mdotF/rf))/(np.exp(-ZT*mdotF/rf)-np.exp(-ZT*mdotF/DelTM))
    ResOuterEnergy = (cpg/LHVheat)*(TGradIn - TGradOut)  -  1
    #Saturation
    A = PsatRef*np.exp(LHVapor/(RFue*TsatRef))
    B = LHVapor/RFue
    ResSaturation = (A*np.exp(-B/Ts)*MWFue)/(A*np.exp(-B/Ts)*MWFue + (Pref-A*np.exp(-B/Ts))*MWPro)  -  YFs
    return [ResInnerMass,ResOuterMass,ResInnerEnergy,ResOuterEnergy,ResSaturation]

def EvapCalc(YOxi,YFue,ReacMech,Pref,Tinf,Tdrop,YOxiInf,rd,rhol,cpl,LHVapor,LHVheat,FAst,TsatRef,PsatRef,MachLiner,Gamma,fracMassEvap):
    gas = ct.Solution(ReacMech)
    #Compute Gas Constant
    Ru = 8.3145 #J/mol/K
    gas.TPY = 298,101325,YFue
    MWFue = gas.mean_molecular_weight/1000 #kg/mole
    RFue =  Ru/MWFue #Gas constant of fuel J/kg/K
    gas.set_equivalence_ratio(1.0,fuel=YFue,oxidizer=YOxi,basis='mass')
    gas.TP = 298,101325
    gas.equilibrate('HP')
    MWPro = gas.mean_molecular_weight/1000
    Tadi = gas.T
    #Initialization
    rdBase = 10e-6 #m radius of droplet
    TinfBase = 487.233
    PrefBase = 485501
    DelTM, dummy1, dummy2, dummy3, dummy4 = DelTMCalc(MachLiner,Gamma,Tinf,Pref,YOxi,Ru,rdBase,gas)
    InitGuess = [1e-9,rdBase*1.1,Tadi,Tdrop,0.5]
    extArgs = [YOxi,YFue,gas,Pref,Tinf,Tdrop,YOxiInf,rd,DelTM,cpl,LHVapor,LHVheat,FAst,TsatRef,PsatRef,MWFue,RFue,MWPro]
    while abs(rdBase-rd)/rdBase >0.0001:
        extArgs[7] = rdBase
        extArgs[8], dummy1, dummy2, dummy3, dummy4 = DelTMCalc(MachLiner,Gamma,TinfBase,PrefBase,YOxi,Ru,rdBase,gas)
        extArgs[3] = PrefBase
        extArgs[4] = TinfBase
        Result = fsolve(Objective,InitGuess,args = extArgs)
        rdBase = (rd-rdBase)/100 + rdBase
        TinfBase = (Tinf-TinfBase)/100 + TinfBase
        PrefBase = (Pref-PrefBase)/100 + PrefBase
        InitGuess = Result
    #Droplet Lifetime Integrator
    fracrdEvap = fracMassEvap**(1/3) #Stop Criterion for radius
    time = np.linspace(0,1e-2,100000)
    rdLst = np.ones(len(time))*rd
    DelTMLst = np.ones(len(time))
    NuLst = np.ones(len(time))
    ReLst = np.ones(len(time))
    PrLst = np.ones(len(time))
    velRelLst = np.ones(len(time))
    mdotFLst = np.ones(len(time))
    rFLst = np.ones(len(time))
    TFLst = np.ones(len(time))
    TSLst = np.ones(len(time))
    YSLst = np.ones(len(time))
    for i in range(1,len(time)):
        rdCur = rdLst[i-1]
        extArgs[7] = rdCur
        DelTMLst[i-1],NuLst[i-1],ReLst[i-1],PrLst[i-1],velRelLst[i-1] = DelTMCalc(MachLiner,Gamma,Tinf,Pref,YOxi,Ru,rdCur,gas)
        extArgs[8] = DelTMLst[i-1]
        extArgs[3] = Pref
        extArgs[4] = Tinf
        Result = fsolve(Objective,InitGuess,args = extArgs)
        InitGuess = Result
        mdotFCur = Result[0]
        #Save History of Vaporization
        mdotFLst[i-1] = Result[0]
        rFLst[i-1] = Result[1]
        TFLst[i-1] = Result[2]
        TSLst[i-1] = Result[3]
        YSLst[i-1] = Result[4]
        #Save History of Vaporization
        rdLst[i] = rdCur + ((-mdotFCur)/(4*np.pi*rhol*(rdCur**2)))*(time[i]-time[i-1])
        if rdLst[i] < rd*fracrdEvap:
            rdLst = rdLst[:(i+1)]
            time = time[:(i+1)]
            DelTMLst = DelTMLst[:(i+1)]
            NuLst = NuLst[:(i+1)]
            ReLst = ReLst[:(i+1)]
            PrLst = PrLst[:(i+1)]
            velRelLst = velRelLst[:(i+1)]
            mdotFLst = mdotFLst[:(i+1)]
            rFLst = rFLst[:(i+1)]
            TFLst = TFLst[:(i+1)]
            TSLst = TSLst[:(i+1)]
            YSLst = YSLst[:(i+1)]
            print("Time Found")
            break
    return [rdLst,time,DelTMLst,NuLst,ReLst,PrLst,velRelLst,mdotFLst,rFLst,TFLst,TSLst,YSLst]