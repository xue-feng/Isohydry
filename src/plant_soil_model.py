
''' Author: Xue Feng; Date: 02/01/2017 '''

import numpy as np
import scipy.optimize as opt
import params_constants 
from params_soil import soil_dict
import scipy.signal as signal
import matplotlib.pyplot as plt

class Canopy:
    def __init__(self, A_canopy, Gs_leaf, c_leaf, Amax): 
        self.A_canopy = A_canopy                    # canopy area
        self.Gs_leaf = Gs_leaf                  # max Gs as function of VPD, e.g. a in a*np.exp(-b*x)
        self.c_leaf = c_leaf                        # parameter for exponential decay as function of P_leaf
        self.Amax = Amax*A_canopy                    # in mol/m2/d, need to convert 
#         self.Amax = Amax*A_canopy/(36000.0*40.5) # need to convert from mol/m2/d to m3/s, per plant
        
    def gmax_canopy(self):
#         gmax_leaf = self.gmax_leaf(VPD)
        gmax_leaf = self.Gs_leaf
        Sd, rhoH2O, Mw = params_constants.Sd, params_constants.rhoH2O, params_constants.Mw
        return gmax_leaf*self.A_canopy*Mw*Sd/(rhoH2O)
    
    def g_canopy_exponential(self, P_leaf):
        # gives canopy conductance g_canopy based on linear decay function 
        gmax = self.gmax_canopy()
        return gmax*np.exp(self.c_leaf*P_leaf)
    
    def g_canopy(self, P_leaf):
        # returns canopy conductance based on specified functional form
        g = self.g_canopy_exponential(P_leaf)
        return g
    
    def carboxyl_constant(self):
        ## uses ca = 400 ppm = 400E-6 mol/mol
        gmax, Amax = self.gmax_canopy(), self.Amax/(36000.0*40.5)
        return -Amax/(Amax*1.6 - gmax*0.000400)
    
    def A(self, P_leaf):
        gs, k = self.g_canopy(P_leaf), self.carboxyl_constant()
        A_inst = gs*k*0.000400/(1.6*k+gs)
        return A_inst*(40.5*36000.0)  #daily assimilation per plant
    
class Stem: 
    def __init__(self, L_stem, A_stem, Ksat_stem, P50_stem, a_stem, plc_form=None):
        self.L_stem = L_stem            # length of conducting stem
        self.A_stem = A_stem            # sapwood area 
        self.Ksat_stem = Ksat_stem      
        self.a_stem = a_stem
        self.P50_stem = P50_stem
        self.plc_form = plc_form
        
    def gmax_stem(self):
        Sd, rhoH2O = params_constants.Sd, params_constants.rhoH2O
        return self.Ksat_stem*Sd*self.A_stem/(self.L_stem*rhoH2O)
    
    def g_stem_plc(self, P_stem):
        # gives stem conductivity g_stem based on exponential-sigmoidal function 
        gmax = self.gmax_stem()
        return gmax*(1.0-1.0/(1.0+np.exp(self.a_stem*(P_stem - self.P50_stem) )))
    
    def g_stem(self, P_stem):
        g = self.g_stem_plc(P_stem)
        return g

    def Px_solvent(self, P_stem, percent=0.90):
        # helper function for solving for P_stem corresponding to designated percent loss of conductivity
        return self.g_stem(P_stem) - (1.0-percent)*self.gmax_stem()
    
class Soil_root:
    def __init__(self, L_root, A_root, d_root, soil_type):
        # select soil types: "sand", "loamy_sand", "sandy_loam", "loam", "clay" to set hydrological properties 
        self.soil_type = soil_type 
        self.L_root = L_root        # mean rooting depth, Zr
        self.A_root = A_root        # root area index
        self.d_root = d_root        # fine root diameter
        
    def gmax_root(self):
        Ksat_soil = soil_dict[self.soil_type]['Ksat']
        nu, rhoH2O, g = params_constants.nu, params_constants.rhoH2O, params_constants.g
        return nu*Ksat_soil/(rhoH2O*g)*np.sqrt(self.A_root/(self.d_root*self.L_root))
    
    def g_soil(self, P_soil):
        # gives soil-root conductivity based on soil water potential 
        b = soil_dict[self.soil_type]['b']
        Psat_soil = soil_dict[self.soil_type]['Ps']
        c = 2.0*b+3.0; d =4.0
        return self.gmax_root()*(Psat_soil/P_soil)**((c-d)/b)
    
class Soil:
    def __init__(self, soil_type):
        # select soil types: "sand", "loamy_sand", "sandy_loam", "loam", "clay" to set hydrological properties 
        self.soil_type = soil_type 
        self.Ksat_soil = soil_dict[self.soil_type]['Ksat']
        self.Psat_soil = soil_dict[self.soil_type]['Ps']
        self.b = soil_dict[self.soil_type]['b']
        self.n = soil_dict[self.soil_type]['n']
        
    def s_solver(self, P_soil):
        # returns soil moisture based on soil water potential 
        return (P_soil/self.Psat_soil)**(-1.0/self.b)
    
    def P_soil_solver(self, s):
        # return soil water potential based on soil moisture
        return self.Psat_soil*s**(-self.b)
    
class Whole_plant:
    def __init__(self, species):
        # set name of species and vapor pressure deficit
        self.species = species
    
    @property
    def canopy(self):
        return self._canopy
    @canopy.setter
    def canopy(self, canopy_object):
        self._canopy = canopy_object
        
    @property
    def stem(self):
        return self._stem
    @stem.setter
    def stem(self, stem_object):
        self._stem =  stem_object
        
    @property
    def soil_root(self):
        return self._soil_root
    @soil_root.setter
    def soil_root(self, soil_root_object):
        self._soil_root =  soil_root_object
        
    @property
    def soil(self):
        return self._soil
    @soil.setter
    def soil(self, soil_object):
        self._soil =  soil_object
    
    
    def flux_solvent(self, (ET, P_stem, P_leaf), P_soil, VPD):
        # helper function for solving coupled hydrodynamic system 
        g_soil, g_stem, g_canopy = self.soil_root.g_soil, self.stem.g_stem, self.canopy.g_canopy
        P_baro = params_constants.P_baro
        return (g_canopy(P_leaf)*(VPD/P_baro) - ET, \
                g_stem(P_stem)*(P_stem-P_leaf) - ET, \
                g_soil(P_soil)*g_stem(P_stem)/(g_soil(P_soil)+g_stem(P_stem)) * (P_soil-P_stem) - ET ) 
#                 g_soil(P_soil)*(P_soil-P_stem) - ET )        
    
    def flux_solver(self, P_soil, VPD, init=None):
        if init is None:
            init = (0.000001,P_soil,P_soil-0.1)
        # calculate system of 3 equations describing internal water states (ET, P_stem, P_leaf) 
        # based on an input of soil water potential and VPD
        ET, P_stem, P_leaf = opt.fsolve(self.flux_solvent, init, (P_soil, VPD), xtol=10e-10)
        return ET, P_stem, P_leaf
    
    def flux_solver_s(self, s, VPD, init=None):
        # calculate system of 3 equations describing internal water states (ET, P_stem, P_leaf) 
        # based on an input of soil moisture and VPD
        P_soil = self.soil.P_soil_solver(s)
        ET, P_stem, P_leaf = self.flux_solver(P_soil, VPD, init)
        return ET, P_stem, P_leaf
    
    def get_sCrit(self, plc=0.80):
        Ploss = opt.fsolve(self.stem.Px_solvent, -0.1, args=(plc,))[0]
        return self.soil.s_solver(Ploss)
    
    def get_Pg12(self):
        c_leaf = self.canopy.c_leaf
        Pg12 = 1.0/c_leaf * np.log(0.12)
        return Pg12
    
    def get_fluxes(self, VPD, s):
        s = s[::-1]
        Flux = np.zeros((len(s), 3))
        P_soil = np.zeros(len(s))
        Flux_init = (0.000001,self.soil.P_soil_solver(s[0]),self.soil.P_soil_solver(s[0])-0.1)
        for i, si in enumerate(s): 
            Flux[i,:] = self.flux_solver_s(si, VPD, Flux_init)
            P_soil[i] = self.soil.P_soil_solver(si)  
#         return Flux, P_soil
        return Flux[::-1, :], P_soil[::-1]
    
    def get_fluxes_scalar(self, VPD, s):
        P_soil = self.soil.P_soil_solver(s) 
        Flux_init = (0.000001,P_soil,P_soil-0.1)
        Flux = self.flux_solver_s(s, VPD, Flux_init)
        return Flux, P_soil
    
    def get_Es_params(self, VPD, s):
        Flux, _ = self.get_fluxes(VPD, s) 
        E = Flux[:,0]
        # deriving related parameters; thresholds for linearizing E-s relationship
        Emax_th = 0.95; sst_th = 0.85; sw_th = 0.01
        Emax = Emax_th*max(E)
        sst = s[np.where(E>sst_th*Emax)[0][0]]
        sw = s[np.where(E>sw_th*Emax)[0][0]]
        return sw, sst, Emax
    
    def get_derived_params(self, VPD, s, alpha, n, Ks, sfc, plc=0.80): 
        Es_params = self.get_Es_params(VPD, s)
        Zr = self.soil_root.L_root
        
        ''' to renormalize to per plant basis, soil water storage need to be modified! - matters both for loss and for rain pulse input 
        e.g., instead of Zr*unit ground area, it's Zr*Ar and alpha*Ar
        '''
        Ar = self.soil_root.A_root
        sw = Es_params[0]; sst = Es_params[1]; ETmax = Es_params[2]
#         gam = n*Zr/alpha; eta = ETmax/(n*Zr); k = Ks/(n*Zr)
#         gam = n*Zr*Ar/alpha; eta = ETmax/(n*Zr*Ar); k = Ks/(n*Zr*Ar)
        gam = (n*Zr*Ar)/(alpha*Ar); eta = ETmax/(n*Zr*Ar); k = Ks/(n*Zr*Ar)
        
        sCrit = self.get_sCrit(plc)
        return gam, eta, k, sw, sst, sCrit
    
    def get_sigma(self, VPD, s):
        Flux, P_soil = self.get_fluxes(VPD,s); P_stem = Flux[:,1]
        # find curvature and slope
        diffP = np.ma.masked_greater(P_soil - P_stem, 10) # differential should not be greater than this -- cavitation limits
        # to make diffP even spaced 
#         dP_p1 = diffP[:-1] - diffP[1:]; dP_p2 = dP_p1[:-1] - dP_p1[1:]
        dP_p1 = np.diff(diffP); dP_p2 = np.diff(dP_p1)
        kappa = dP_p2/(1+dP_p1[:-1]**2)**(3.0/2.0) # for curvature
        imax = signal.argrelmax(diffP)[0][-1]
        try: 
            imin = signal.argrelmax(kappa[:-1]-kappa[1:])[0][-4] # going from high moisture, the first minimum curvature
        except IndexError: 
            imin = 0
        sigma = 1.0 - (diffP[imax] - diffP[imin])/(P_soil[imax]-P_soil[imin])
        plt.plot(P_soil, P_stem); plt.plot(P_soil, P_soil, '--'); plt.plot([P_soil[imin],P_soil[imax]],[P_stem[imin],P_stem[imax]], 'o')
#         plt.plot(P_soil, diffP); plt.plot([P_soil[imin], P_soil[imax]], [diffP[imin], diffP[imax]], 'o'); plt.show()
        return sigma
    