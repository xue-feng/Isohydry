''' Author: Xue Feng; Date: 02/01/2017 '''

import numpy as np
import matplotlib.pyplot as plt
import params_constants
from params_soil import soil_dict
import math
import scipy.signal as signal
from utility_functions import import_traits_data, initialize_plant, get_fluxes


def plot_supply_demand(sp, params, s0,  VPD, newfig=True):
    if type(sp) is str:
        plant = initialize_plant(sp, params)
    else: plant = sp
    
    P_soil = plant.soil.P_soil_solver(s0)
    g_soil, g_stem, g_canopy = plant.soil_root.g_soil, plant.stem.g_stem, plant.canopy.g_canopy
    
    Px_arr = np.linspace(-10.0, -0.01, 100)
    P_baro = params_constants.P_baro
    
    ETsupply = g_soil(P_soil)*g_stem(Px_arr)/(g_soil(P_soil)+g_stem(Px_arr)) * (P_soil - Px_arr) 
    Pl_arr = Px_arr - ETsupply/g_stem(Px_arr)
    ETdemand = g_canopy(Pl_arr, VPD)*(VPD/P_baro)
    
    if newfig: plt.figure(figsize=(6,4.5))
    else: plt.figure(3,figsize=(6,4.5))
    plt.plot(Px_arr, ETsupply)
    plt.plot(Px_arr, ETdemand)
    plt.hlines(0, -10.0, 0, linestyle='--')
    


def plot_Es_phaseplots(sp, params, VPD, newfig=True):
    # Returns parameters derived from plant traits: sw, sst, Emax, s_P90, s_Pg12, L_root
    # the argument sp accepts species name as a string, or as a plant class
    if type(sp) is str:
        plant = initialize_plant(sp, params)
    else: plant = sp
    Flux, P_soil = get_fluxes(plant, VPD, s) 
    E = Flux[:,0]; P_stem = Flux[:,1]

    # deriving indices to calculate sigma
    diffP = P_soil - P_stem
    diffmax = np.max(diffP)
    diffmin = 0.01*diffmax
    imax = np.argmax(diffP)
    imin = np.where(diffP<diffmin)[0]
    ilastmin = np.argmax(imin[:-1] - imin[1:] != -1)
    
    if newfig: plt.figure(figsize=(6,4.5))
    else:plt.figure(1,figsize=(6,4.5))
    plt.plot(s, E); plt.ylabel('m/s')
    plt.title('ET vs. soil moisture')
    
    if newfig: plt.figure(figsize=(6,4.5))
    else:plt.figure(2,figsize=(6,4.5))
    plt.plot(s, P_soil-P_stem); plt.ylabel('MPa')
    plt.plot(s[imax], diffP[imax], 'ro')
    plt.plot(s[ilastmin], diffP[ilastmin], 'bo')
    plt.title('Predawn - midday vs. soil moisture')
    
if __name__ =='__main__':
    # import trait parameters
    traits_path = '../traits_data.xls'
    params = import_traits_data() # import chaparral species params
    
    # designate soil type
    soil_type = 'sandy_loam'; smin = soil_dict[soil_type]['sh']
    s = np.arange(smin,1.00,0.001)   
    sp = 'JUNI'; s0=0.156; VPD=1.0; P_baro = params_constants.P_baro
    
    plant = initialize_plant(sp, params)
    P50 = plant.stem.P50_stem
    c_leaf = plant.canopy.c_leaf; Pg12 = 1.0/c_leaf * np.log(0.12)

    shift_range = np.linspace(0.5, 1.5, 20)
    P50_arr = shift_range*P50
    Pg12_arr = shift_range*Pg12
    sigma_grid = np.zeros((len(P50_arr), len(Pg12_arr)))
    for i, pg12 in enumerate(Pg12_arr):
        for j, p50 in enumerate(P50_arr):
            print i, j, pg12, p50 
            par_modified = params.copy()
            par_modified[sp]['stem_dict']['P50_stem'] = p50
            par_modified[sp]['canopy_dict']['c_leaf'] = c_leaf/(pg12/Pg12)
            plant = initialize_plant(sp, par_modified)
            Flux = np.zeros((len(s), 3))
            P_soil = np.zeros(len(s))
            for k, si in enumerate(s):
                P_soil[k] = plant.soil.P_soil_solver(si) 
                Flux[k,:] = plant.flux_solver_s(si, VPD) 
            P_stem = Flux[:,1]; P_leaf = Flux[:,2]
            diffP = np.ma.masked_greater(P_soil - P_stem, 10)
            diffmax = np.max(diffP)
            
            dP_p1 = diffP[:-1] - diffP[1:]
            dP_p2 = dP_p1[:-1] - dP_p1[1:]
            kappa = dP_p2/(1+dP_p1[:-1]**2)**(3.0/2.0)
            imax = signal.argrelmax(diffP)[0][-1]
            try: 
                imin = signal.argrelmax(kappa[:-1]-kappa[1:])[0][-2] # going from high moisture, the first minimum curvature
                sigma_grid[i,j] = math.atan(np.pi/4.0 - np.tan(diffP[imax]/(P_soil[imax]-P_soil[imin])))
            except IndexError:  
                sigma_grid[i,j] = np.nan
    PX50, PL12 = np.meshgrid(P50_arr, Pg12_arr)
    plt.figure(figsize=(4.8,4.5))
    CS = plt.contour(PX50, PL12, sigma_grid)
    plt.clabel(CS, fontsize=10, inline=1)
    plt.show()
    
    
    plant = initialize_plant(sp, params)
    P_soil = plant.soil.P_soil_solver(s0)
    g_stem, g_soil, g_canopy = plant.stem.g_stem, plant.soil_root.g_soil, plant.canopy.g_canopy
    print g_soil(P_soil)
    
    Px_arr = np.linspace(P_soil-10.0, P_soil, 100)
    ETsupply = g_soil(P_soil)*g_stem(Px_arr)/(g_soil(P_soil)+g_stem(Px_arr)) * (P_soil - Px_arr) 
    Pl_arr = Px_arr - ETsupply/g_stem(Px_arr)
    ETdemand = g_canopy(Pl_arr, VPD)*(VPD/P_baro)
    
#     plt.plot(Px_arr, g_canopy(Px_arr, VPD))
#     plt.plot(Px_arr, Pl_arr)
#     plt.plot(plant.soil.P_soil_solver(s), g_soil(plant.soil.P_soil_solver(s)))
#     plt.plot(Px_arr, ETsupply)
#     plt.plot(Px_arr, ETdemand)
#     plt.plot(Px_arr, g_stem(Px_arr), '--')
#     plt.plot(Px_arr,g_soil(P_soil)*g_stem(Px_arr)/(g_soil(P_soil)+g_stem(Px_arr))) 
# plt.plot(plant.soil.P_soil_solver(s), g_soil(plant.soil.P_soil_solver(s)))
# plt.yscale('log')
# plt.xlim(-20,0)
# plt.show()

#     plot_Es_phaseplots(sp, params, VPD, newfig=False)
#     plot_supply_demand(sp, params, s0, VPD, newfig=False)
    
    ''' tweak the parameters a bit to test hypothesis '''
    ''' the stem parameters have much lower effect... why??? '''
#     params[sp]['stem_dict']['A_stem'] = 2*params[sp]['stem_dict']['A_stem']
#     params[sp]['stem_dict']['a_stem'] -= 0.3
#     params[sp]['stem_dict']['P50_stem'] += 4.8
#     params[sp]['root_dict']['L_root'] += 2.0
    params[sp]['canopy_dict']['b_leaf'] -= 0.3  
#     soil_type = 'sandy_loam'; soil = Soil(soil_type)
#     params[sp]['canopy_dict']['c_leaf'] -= 0.3
#     params[sp]['canopy_dict']['A_canopy'] = 0.7*params[sp]['canopy_dict']['A_canopy']
#     params[sp]['stem_dict']['Ksat_stem'] = 2.0 * params[sp]['stem_dict']['Ksat_stem'] 
#     s0 -= 0.01
    
    plant = initialize_plant(sp, params)
    g_stem, g_soil, g_canopy = plant.stem.g_stem, plant.soil_root.g_soil, plant.canopy.g_canopy
    P_soil = plant.soil.P_soil_solver(s0)
    
    ETsupply = g_soil(P_soil)*g_stem(Px_arr)/(g_soil(P_soil)+g_stem(Px_arr)) * (P_soil - Px_arr) 
    Pl_arr = Px_arr - ETsupply/g_stem(Px_arr)
    ETdemand = g_canopy(Pl_arr, VPD)*(VPD/P_baro)
      
#     plt.plot(Px_arr, g_canopy(Px_arr, VPD))
#     plt.plot(Px_arr, Pl_arr)
#     plt.plot(Px_arr, (P_soil - Px_arr))
#     plt.plot(Px_arr, g_stem(Px_arr), '--')
#     plt.plot(Px_arr, g_soil(P_soil)*g_stem(Px_arr)/(g_soil(P_soil)+g_stem(Px_arr))) 
#     plt.plot(Px_arr, ETsupply )
#     plt.plot(Px_arr, ETdemand )
#     plt.plot(plant.soil.P_soil_solver(s), g_soil(plant.soil.P_soil_solver(s))); plt.yscale('log'); plt.xlim(-5,0)

#     plot_Es_phaseplots(sp, params, VPD, newfig=False)
#     plot_supply_demand(sp, params, s0, VPD, newfig=False)
    
    plt.show()
    