''' Author: Xue Feng; Date: 02/01/2017 '''

import numpy as np
import matplotlib.pyplot as plt
import params_constants
from params_soil import soil_dict
import copy
import scipy.signal as signal
import math
from utility_functions import import_traits_data, initialize_plant


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
    Flux, P_soil = get_fluxes(plant, VPD,s) 
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

def plot_Es_trait_dependence(trait_name, part_name, trait_arr, plot_arr):
    par_modified = copy.deepcopy(params)
    EsParams_arr = np.zeros((len(shift_range), 3))
    for i, trait_val in enumerate(trait_arr): 
        print i, trait_val
        par_modified[sp][part_name][trait_name] = trait_val
        plant = initialize_plant(sp, par_modified)
        EsParams_arr[i,:] = get_Es_params(plant,VPD,s)  
#     plt.plot(EsParams_arr[:,0], EsParams_arr[:,1], 'o'); plt.xlim(0.15, 0.22); plt.ylim(0.25,0.45); plt.show()
#     plt.figure(figsize=(6,4.5))
#     plt.plot(plot_arr, EsParams_arr[:,0]) # sw
#     plt.plot(plot_arr, EsParams_arr[:,1]) # sst
#     plt.plot(plot_arr, EsParams_arr[:,2]) # Emax
#     plt.plot(plot_arr, EsParams_arr[:,1] - EsParams_arr[:,0])
#     plt.ylim(0,sfc)

def plot_sigma_trait_dependence(trait_name, part_name, trait_arr, plot_arr):
    par_modified = copy.deepcopy(params)
    sigma_arr = np.zeros(len(shift_range))
    for i, trait_val in enumerate(trait_arr): 
        print i, trait_val,
        par_modified[sp][part_name][trait_name] = trait_val
        plant = initialize_plant(sp, par_modified)
        Flux, P_soil = get_fluxes(plant, VPD,s); P_stem = Flux[:,1]
        # find curvature and slope
        diffP = np.ma.masked_greater(P_soil - P_stem, 10) # differential should not be greater than this -- cavitation limits
        # to make diffP even spaced 
        dP_p1 = diffP[:-1] - diffP[1:]
        dP_p2 = dP_p1[:-1] - dP_p1[1:]
        kappa = dP_p2/(1+dP_p1[:-1]**2)**(3.0/2.0)
        imax = signal.argrelmax(diffP)[0][-1]
        try: 
            imin = signal.argrelmax(kappa[:-1]-kappa[1:])[0][-4] # going from high moisture, the first minimum curvature
            print imin
            sigma_arr[i] = math.atan(np.pi/4.0 - np.tan(diffP[imax]/(P_soil[imax]-P_soil[imin])))
#             plt.plot(P_soil, diffP); plt.plot([P_soil[imin], P_soil[imax]], [diffP[imin], diffP[imax]], 'o')
        except IndexError: 
            sigma_arr[i] = np.nan
        
#     plt.show() 
#     plt.figure(figsize=(6,4.5))
    plt.plot(plot_arr, sigma_arr)
    plt.ylim(0.5,0.7)
    
if __name__ =='__main__':
    # import trait parameters
    traits_path = '../traits_data.xls'
    params = import_traits_data() # import chaparral species params
    
    # designate soil type
    soil_type = 'sandy_loam';  smin = soil_dict[soil_type]['sh']
    sfc = soil_dict[soil_type]['sfc']; sst = soil_dict[soil_type]['sst']
#     s = np.linspace(smin, 1.0, 1000)
    s = np.hstack((np.linspace(smin,sst,1000), np.linspace(sst+0.001,1.0,100)))
    sp = 'PINE'; s0=0.156; VPD=2.0; P_baro = params_constants.P_baro
    print get_Es_params(initialize_plant(sp, params), VPD,s)
    
    # setting up parameter ranges 
    plant = initialize_plant(sp, params); shift_range = np.linspace(0.1, 2.0, 30)
    
    P50 = plant.stem.P50_stem
    c_leaf = plant.canopy.c_leaf; Pg12 = 1.0/c_leaf * np.log(0.12)
    
    A_canopy = plant.canopy.A_canopy
    Zr = plant.soil_root.L_root
    P50_arr = shift_range*P50
    
    Pg12_arr = shift_range*Pg12
    Acan_arr = shift_range*A_canopy
    Zr_arr = shift_range*Zr
    
#     plt.figure(figsize=(6,4.5))
#     plot_Es_trait_dependence('P50_stem', 'stem_dict', P50_arr, P50_arr)
#     plot_Es_trait_dependence('c_leaf', 'canopy_dict', c_leaf/shift_range, Pg12_arr)
#     plot_Es_trait_dependence('A_canopy', 'canopy_dict', Acan_arr, Acan_arr)
#     plot_Es_trait_dependence('L_root', 'root_dict', Zr_arr, Zr_arr)
#     plt.ylim(0.10, 0.25)
#     plt.show()
    
    plt.figure(figsize=(6,4.5))
    plot_sigma_trait_dependence('P50_stem', 'stem_dict', P50_arr, P50_arr)
    plot_sigma_trait_dependence('c_leaf', 'canopy_dict', c_leaf/shift_range, Pg12_arr)
    plot_sigma_trait_dependence('A_canopy', 'canopy_dict', Acan_arr, Acan_arr)
    plot_sigma_trait_dependence('L_root', 'root_dict', Zr_arr, Zr_arr)
    plt.ylim(0.55,0.70)
    plt.show()       