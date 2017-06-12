
''' check whether simulated parameters produce behaviors were realistic according to Mencuccini '''

import matplotlib.pyplot as plt
import cPickle as pickle
import params_constants
import numpy as np
from scipy import stats
from params_soil import soil_dict
from SALib.analyze import sobol
from utility_functions import import_traits_data, get_part, initialize_generic_plant, simulate_ps_t
import matplotlib.ticker as mtick

## soil conditions ##
soil_type = 'loamy_sand'; soil_params = soil_dict[soil_type]
smin = soil_params['sh']
sfc = soil_params['sfc']
sst = soil_params['sst']
n = soil_params['n']
Ks = soil_params['Ksat']
Ps = soil_params['Ps']
b = soil_params['b']
traits = import_traits_data()
s_sat = (-0.005/Ps)**(-1.0/b)
s = np.hstack((np.linspace(smin,sst,500), np.linspace(sst+0.001,1.0,100)))

n_trajectories = 100; dt = 0.1
lam=0.15; alpha=0.010; s1 = sfc; s0 = 0.5; Amax = 1.0/dt; Rmax = 0.10*Amax; plc = 0.5

## prep for sensitivity analysis ##
var_names = np.array(['A_canopy','Gs_leaf','c_leaf','L_stem','A_stem','Ksat_stem','a_stem','P50_stem','L_root','A_root','Rmax'])
n_vars = len(var_names)

sp = 'JUNI'; tmax = 180; VPD=2.0; tRun = np.arange(0,tmax+dt,dt)

with open('../Si_'+sp+'_vpd'+str(int(VPD))+'_tmax'+str(tmax)+'_severe_params.pickle', 'rb') as handle:
    params = pickle.load(handle)
count1 = []
count2 = []
print 'total number of params:', len(params)
for i, p_vals in enumerate(params):
    print i
    plant_traits = p_vals[:-1]; R = p_vals[-1]
    generic_plant = initialize_generic_plant(var_names, plant_traits, soil_type)
#     (ET, P_stem, P_leaf), Psoil = generic_plant.get_fluxes_scalar(VPD, s_sat)
#     if (P_leaf<-3.0) :  ## Criteria 1
#         count1 += 1
#         print i, P_leaf
    
    # Check criteria 2: mean soil moisture, assimilation/stomatal conductance did not decrease during dry down
    gam,eta,k,sw,sst,sCrit=generic_plant.get_derived_params(VPD, s, alpha, n, Ks, sfc, plc=plc)
    ps, assm = simulate_ps_t(n_trajectories, tRun, dt, s0, lam, gam, eta, k, sw, sst, s1, Amax, R)
    
    assm_t = np.mean(assm+R,axis=0)
    s_t = np.mean(ps,axis=0)
    gs = np.zeros(len(s_t))
    for j, si in enumerate(s_t): 
        (ET, P_stem, P_leaf), Psoil = generic_plant.get_fluxes_scalar(VPD, si)
        gs[j] = generic_plant.canopy.g_canopy(P_leaf, VPD)
        
    gs_re = np.reshape(gs[1:], (tmax,10))
    assm_re = np.reshape(assm_t[1:], (tmax,10))
    assm_gs_ratio = np.mean(assm_re, axis=1)/np.mean(gs_re,axis=1)
    if (np.argmax(assm_gs_ratio) < 0.9*len(assm_gs_ratio)) and (assm_gs_ratio[-1]<0.8*(np.max(assm_gs_ratio)-np.min(assm_gs_ratio))):
        count2.append(i)
        print i, np.argmax(assm_gs_ratio)
        plt.plot(assm_gs_ratio); plt.plot(np.mean(assm_re, axis=1)); plt.plot(np.mean(gs_re, axis=1))

# print 'criteria 1 failed', float(count1)/len(params)

print 'criteria 2 failed', float(len(count2))/len(params)   
print count2
plt.show()