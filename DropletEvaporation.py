import numpy as np
import cantera as ct
from scipy.optimize import fsolve
#DesignParameters
Tinf = 1000#K
Pref = 201325#Pa
Yvinf = 0#
rd = 100*1e-6 #50 micron
rinf = rd*10 #1cm
numPt = 100
LHVv = 335*1e3 #J/kg
hfg = 256e3 #J/kg
Psatref = 101325 #Pa
Tsatref = 489.45 #K
#Initial Conditions
YAir = 'O2:0.2329,N2:0.7671'
YFuel = 'H2O:1'
gas = ct.Solution('gri30.yaml')
gas.TPY = 298,101325,YAir
YAir = gas.Y
gas.TPY = 298,101325,YFuel
YFuel = gas.Y
rLst = np.linspace(rd,rinf,numPt)
TLstIni = (rLst-rLst[0])/(rLst[-1]-rLst[0])*(Tinf*0.5) + Tinf*0.5
YvLstIni = 1 - (rLst-rLst[0])/(rLst[-1]-rLst[0])*(1-Yvinf)
mdotIni = 1e-15
InitGuess = np.concatenate((TLstIni,YvLstIni,np.array([mdotIni])))
#Objective Function
def ThermalUpdator(gas,YvLst,TLst,Pref,YFuel,YAir):
    DiffLst = np.zeros(len(TLst))
    rhoLst = np.zeros(len(TLst))
    condLst = np.zeros(len(TLst))
    hLst = np.zeros(len(TLst))
    for idx in range(len(TLst)):
        gas.TPY = TLst[idx],Pref,(YFuel*YvLst[idx]+YAir*(1-YvLst[idx]))
        DiffLst[idx] = gas.binary_diff_coeffs[np.argmax(YFuel)][np.argmax(YAir)]
        rhoLst[idx] = gas.density_mass
        condLst[idx] = gas.thermal_conductivity
        hLst[idx] = gas.enthalpy_mass
    return DiffLst,rhoLst,condLst,hLst

def PsatCalc(gas,Ti,Pref,YFuel,YAir,Psatref,Tsatref,hfg):
    Ru = 8.3145 #J/mol/K
    gas.TPY = Ti,Pref,YFuel
    MWFuel = gas.mean_molecular_weight/1000 #kg/mole
    gas.TPY = Ti,Pref,YAir
    MWAir = gas.mean_molecular_weight/1000 #kg/mole
    R = Ru/MWFuel #J/kg/K
    C1 = np.log(Psatref)+(hfg/R)*(1/Tsatref)
    Psat = np.exp(C1-(hfg/(R*Ti)))
    Ysat = 1/(1+(MWAir/MWFuel)*((Pref/Psat)-1))
    return Ysat

def Objective(vars,extArgs):
    rLst  = extArgs[0]
    YFuel = extArgs[1]
    YAir  = extArgs[2]
    Tinf  = extArgs[3]
    Pref  = extArgs[4]
    Yvinf = extArgs[5]
    gas   = extArgs[6]
    Psatref = extArgs[7]
    Tsatref = extArgs[8]
    hfg = extArgs[9]
    TLst  = vars[:len(rLst)]
    YvLst = vars[len(rLst):-1]
    mdot  = vars[-1]
    #Update thermal variables
    DiffLst,rhoLst,condLst,hLst = ThermalUpdator(gas,YvLst,TLst,Pref,YFuel,YAir)
    rhoDLst = DiffLst*rhoLst
    rhoD_ip1 = rhoDLst[3:]
    rhoD_i   = rhoDLst[2:-1]
    rhoD_im1 = rhoDLst[1:-2]
    cond_ip1 = condLst[3:]
    cond_i   = condLst[2:-1]
    cond_im1 = condLst[1:-2]
    velLst = mdot/(4*np.pi*(rLst**2)*rhoLst)
    htLst = hLst+0.5*(velLst**2)
    ht_i   = htLst[2:-1]
    ht_im1 = htLst[1:-2]
    ht_im2 = htLst[0:-3]
    #Species Equation
    Delr = np.diff(rLst)
    Yv_ip1 = YvLst[3:]
    Yv_i   = YvLst[2:-1]
    Yv_im1 = YvLst[1:-2]
    Yv_im2 = YvLst[0:-3]
    Delr_i   = Delr[2:]
    Delr_im1 = Delr[1:-1]
    Delr_im2 = Delr[0:-2]
    YadvMain = mdot*((3*Yv_i-4*Yv_im1+Yv_im2)/(Delr_im2+Delr_im1))
    r_i      = rLst[2:-1]
    d2Yvdr2Main = (Yv_ip1-Yv_i)/(Delr_i**2) - (Yv_i-Yv_im1)/(Delr_i*Delr_im1)
    dYvdrMain   = (Yv_ip1-Yv_im1)/(Delr_im1+Delr_i)
    drhoDdrMain = (rhoD_ip1-rhoD_im1)/(Delr_im1+Delr_i)
    YdifMain    = -4*np.pi*(rhoD_i*(r_i**2)*d2Yvdr2Main + 2*r_i*rhoD_i*dYvdrMain + (r_i**2)*dYvdrMain*drhoDdrMain)
    ResYMain    = YadvMain + YdifMain
    #Energy Equation
    T_ip1 = TLst[3:]
    T_i   = TLst[2:-1]
    T_im1 = TLst[1:-2]
    TadvMain    = mdot*((3*ht_i-4*ht_im1+ht_im2)/(Delr_im2+Delr_im1))
    d2Tdr2Main  = (T_ip1-T_i)/(Delr_i**2) - (T_i-T_im1)/(Delr_i*Delr_im1)
    dTdrMain    = (T_ip1-T_im1)/(Delr_im1+Delr_i)
    dconddrMain = (cond_ip1-cond_im1)/(Delr_im1+Delr_i)
    TdifMain    = -4*np.pi*(cond_i*(r_i**2)*d2Tdr2Main + 2*r_i*cond_i*dTdrMain + (r_i**2)*dTdrMain*dconddrMain)
    ResTMain    = TadvMain + TdifMain
    #Special treatment of node 1 Y (to the immediate right of droplet surface)
    YadvSpe    = mdot*((YvLst[1]-YvLst[0])/(Delr[0]))
    d2Yvdr2Spe = (YvLst[2]-YvLst[1])/(Delr[1]**2) - (YvLst[1]-YvLst[0])/(Delr[1]*Delr[0])
    dYvdrSpe   = (YvLst[2]-YvLst[0])/(Delr[0]+Delr[1])
    drhoDdrSpe = (rhoDLst[2]-rhoDLst[0])/(Delr[0]+Delr[1])
    YdifSpe    = -4*np.pi*(rhoDLst[1]*(rLst[1]**2)*d2Yvdr2Spe + 2*rLst[1]*rhoDLst[1]*dYvdrSpe + (rLst[1]**2)*dYvdrSpe*drhoDdrSpe)
    ResYSpe    = YadvSpe + YdifSpe
    #Special treatment of node 1 T (to the immediate right of droplet surface)
    TadvSpe    = mdot*((htLst[1]-htLst[0])/(Delr[0]))
    d2Tdr2Spe  = (TLst[2]-TLst[1])/(Delr[1]**2) - (TLst[1]-TLst[0])/(Delr[1]*Delr[0])
    dTdrSpe    = (TLst[2]-TLst[0])/(Delr[0]+Delr[1])
    dconddrSpe = (condLst[2]-condLst[0])/(Delr[0]+Delr[1])
    TdifSpe    = -4*np.pi*(condLst[1]*(rLst[1]**2)*d2Tdr2Spe + 2*rLst[1]*condLst[1]*dTdrSpe + (rLst[1]**2)*dTdrSpe*dconddrSpe)
    ResTSpe    = TadvSpe + TdifSpe
    #Condtion 1 Ts->Ys
    ResTs2Ys = YvLst[0] - PsatCalc(gas,TLst[0],Pref,YFuel,YAir,Psatref,Tsatref,hfg)
    #Condition 2 Ts->mdot
    ResTs2Mdot = mdot*hfg - 4*np.pi*(rLst[0]**2)*condLst[0]*((TLst[1]-TLst[0])/Delr[0])
    #Condition 3 Ys->mdot
    ResYs2Mdot = mdot*YvLst[0] - 4*np.pi*(rLst[0]**2)*rhoDLst[0]*((YvLst[1]-YvLst[0])/Delr[0]) - mdot
    #Condition 4&5 Far Field
    ResTsFar = Tinf - TLst[-1]
    ResYsFar = Yvinf - YvLst[-1]
    ResFinal = np.concatenate((ResYMain,ResTMain,np.array([ResYSpe,ResTSpe,ResTs2Ys,ResTs2Mdot,ResYs2Mdot,ResTsFar,ResYsFar])))
    return ResFinal

Result = fsolve(Objective,InitGuess,args = [rLst,YFuel,YAir,Tinf,Pref,Yvinf,gas,Psatref,Tsatref,hfg])[0] #Reference Mach number & try initial guess between [0 and 1]









    



