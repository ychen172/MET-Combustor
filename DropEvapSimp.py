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