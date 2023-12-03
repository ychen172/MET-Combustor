import numpy as np
import cantera as ct
from scipy.optimize import fsolve
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

#DesignParameters
YOxi = "O2:0.2314,N2:0.7622,H2O:0.0064"
YFue = "C2H5OH:1"
gas = ct.Solution('CRECK2019-NOx-NoBin.yaml')
Pref = 101325*15 #Pa
Tinf = 800 #K
Tdrop = 298 #K Droplet Temperature
YOxiInf = 1.0 #Farfield Composition
rd = 25e-6 #m radius of droplet
DelTM = 3 #Unity Lewis Number Nu = Sh and so DelT = DelM
rhol = 789 #kg/m3
cpl = 2570 #J/kg/K
LHVapor = 918187.9354880721 #J/kg
LHVheat = 27728322.542190198 #J/kg
FAst = 0.11107613244759687 #Stoichiometric Fuel Air Ratio
TsatRef = 298 #K Fuel saturation temperature
PsatRef = 0.008e6 #Pa 0.008 MPa 
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
#Solve
InitGuess = [1e-10,rd*1.04,Tadi,298,0.9]
extArgs = [YOxi,YFue,gas,Pref,Tinf,Tdrop,YOxiInf,rd,DelTM,cpl,LHVapor,LHVheat,FAst,TsatRef,PsatRef,MWFue,RFue,MWPro]
Result = fsolve(Objective,InitGuess,args = extArgs)[0] 
