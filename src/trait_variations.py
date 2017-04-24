from utility_functions import import_traits_data, initialize_plant, simulate_ps_t
from params_soil import soil_dict
import numpy as np
import matplotlib.pyplot as plt
import copy

def initialize_modified_plant(sp, soil_type, part_name, trait_name, trait_val):
    traits_modified = copy.deepcopy(traits)
    traits_modified[sp][part_name][trait_name] = trait_val
    plant_modified = initialize_plant(sp, traits_modified, soil_type)
    return plant_modified

def get_psCrit(ps, sCrit):
    return len(ps[ps<sCrit])/float(np.shape(ps)[0]*np.shape(ps)[1])

def get_relGnet(assm, R):
    return np.mean(np.sum(assm, axis=1))/((Amax-R)*np.shape(assm)[1]) 
    
def get_modified_params(sp, part_name, trait_name, trait_arr, param_name, VPD):
    params_arr = np.zeros((len(trait_arr), 3)) # allow for three slots to fill up to three params (e.g., Emax, sw, sst)
    if trait_name == 'R': plant = sp
    for i, trait_val in enumerate(trait_arr): 
        print i, trait_val
        if trait_name != 'R': 
            modified_plant = initialize_modified_plant(sp, soil_type, part_name, trait_name, trait_val)
        if param_name=='Es': 
            func = modified_plant.get_Es_params
            params_arr[i,:] = func(VPD, s)
        elif param_name=='sigma': 
            func = modified_plant.get_sigma
            params_arr[i,:] = func(VPD, s)
        elif param_name=='Cgain' or param_name=='Hrisk':
            if trait_name == 'R':
                modified_plant = plant
                gam,eta,k,sw,sst,sCrit=modified_plant.get_derived_params(VPD, s, alpha, n, Ks, sfc, plc=plc)
                ps, assm = simulate_ps_t(n_trajectories, tRun, dt, s0, lam, gam, eta, k, sw, sst, s1, Amax, trait_val)
            else: 
                gam,eta,k,sw,sst,sCrit=modified_plant.get_derived_params(VPD, s, alpha, n, Ks, sfc, plc=plc)
                ps, assm = simulate_ps_t(n_trajectories, tRun, dt, s0, lam, gam, eta, k, sw, sst, s1, Amax, Rmax)
            params_arr[i,0] = get_psCrit(ps, sCrit)
            params_arr[i,1] = get_relGnet(assm, Rmax)  
    return params_arr
 
def plot_param_trait_dependence(plant, param_name, VPD):
    # setting up parameter ranges 
    shift_range = np.linspace(0.1, 2.0, 10)
    sp = plant.species
    
    P50 = plant.stem.P50_stem
    c_leaf = plant.canopy.c_leaf; Pg12 = 1.0/c_leaf * np.log(0.12)
    A_canopy = plant.canopy.A_canopy
    Zr = plant.soil_root.L_root
    Gs_leaf = plant.canopy.Gs_leaf
    
    P50_arr = shift_range*P50
    Pg12_arr = shift_range*Pg12
    Acan_arr = shift_range*A_canopy
    Zr_arr = shift_range*Zr
    Gs_arr = shift_range*Gs_leaf
    R_arr = shift_range*Rmax
    
    plt.figure(figsize=(6,4.5))
    Params_R = get_modified_params(plant, '', 'R', R_arr, param_name, VPD)
    Params_cleaf = get_modified_params(sp, 'canopy_dict','c_leaf', c_leaf/shift_range, param_name, VPD)
    Params_P50 = get_modified_params(sp,'stem_dict', 'P50_stem', P50_arr, param_name, VPD)
    Params_Acanopy = get_modified_params(sp, 'canopy_dict','A_canopy', Acan_arr, param_name, VPD)
    Params_Lroot = get_modified_params(sp,'root_dict','L_root', Zr_arr, param_name, VPD)
    Params_Gs = get_modified_params(sp,'canopy_dict','Gs_leaf', Gs_arr, param_name, VPD)
    
    if param_name == 'Es': ip=2
    elif param_name == 'sigma': ip=0
    elif param_name == 'Hrisk': ip=0
    elif param_name == 'Cgain': ip=1
    plt.plot(P50_arr, Params_P50[:,ip], lw=1, label='P50') 
    plt.plot(Pg12_arr, Params_cleaf[:,ip], lw=1, label='Pg12') 
    plt.plot(Acan_arr, Params_Acanopy[:,ip], lw=1, label='Acanopy') 
    plt.plot(Zr_arr, Params_Lroot[:,ip], lw=1, label='Zr') 
    plt.plot(Gs_arr, Params_Gs[:,ip], lw=1, label='Gs') 
    plt.plot(R_arr, Params_R[:,ip], lw=1, label='R') 
    plt.legend(loc=1)

    
# set soil conditions 
soil_type = 'loamy_sand'
smin = soil_dict[soil_type]['sh']
sfc = soil_dict[soil_type]['sfc']
sst = soil_dict[soil_type]['sst']
n = soil_dict[soil_type]['n']
Ks = soil_dict[soil_type]['Ksat']
s = np.hstack((np.linspace(smin,sst,500), np.linspace(sst+0.001,1.0,100)))

# set species traits
traits = import_traits_data()
juni_plant = initialize_plant('JUNI', traits, soil_type)
pine_plant = initialize_plant('PINE', traits, soil_type)

n_trajectories = 200; dt = 0.1
lam=0.15; alpha=0.010; s1 = sfc; s0 = 0.5; Amax = 1.0/dt; Rmax = 0.10*Amax; plc = 0.5
VPD = 2.0; tmax = 180; tRun = np.arange(0,tmax+dt,dt)

# plot_param_trait_dependence(juni_plant, 'Hrisk', VPD=2.0)
# # plot_param_trait_dependence(juni_plant, 'sigma', VPD=2.0)
# plt.show()

def plot_metrics_param_dependence(plant, VPD, tmax, params='beta'):
    '''metrics: (both are plotted)
    hrisk: hydraulic risk
    assim: carbon assimilation
    params: 
    alpha (structual/morphological traits), 
    beta (risk aversity traits), 
    kappa (transport efficiency),
    sigma (ratio of differential to predawn)
    
    '''
    # setting up parameter ranges 
    shift_range = np.linspace(0.1, 2.0, 10)
    sp = plant.species
    # getting plant traits
    P50 = plant.stem.P50_stem
    c_leaf = plant.canopy.c_leaf; # Pg12 = 1.0/c_leaf * np.log(0.12)
    A_canopy = plant.canopy.A_canopy
    Zr = plant.soil_root.L_root
    Gs_leaf = plant.canopy.Gs_leaf
    A_root = plant.soil_root.A_root
    Ksat_stem = plant.stem.Ksat_stem
    # constructing range of traits
    P50_arr = shift_range*P50
    Acan_arr = shift_range*A_canopy
    Zr_arr = shift_range*Zr
    Gs_arr = shift_range*Gs_leaf
    R_arr = shift_range*Rmax
    Ksat_arr = shift_range*Ksat_stem
    # set up simulation
    param_ranges = np.array([P50_arr, c_leaf/shift_range, Acan_arr, Zr_arr, Gs_arr, R_arr, Ksat_arr])
    part_names = np.array(['stem_dict', 'canopy_dict', 'canopy_dict', 'root_dict', 'canopy_dict', '', 'stem_dict'])
    trait_names = np.array(['P50_stem', 'c_leaf', 'A_canopy', 'L_root', 'Gs_leaf', 'R', 'Ksat_stem'])
    colors = np.array(['orange', 'blue', 'green', 'brown', 'purple', 'black', 'red'])
    # select subset under investigation 
    if params=='alpha': ip = [3,2]
    elif params=='beta': ip = [0,1]
    elif params=='kappa': ip = [4]
    elif params=='phi': ip=[5]
    elif params=='sigma': ip = [0,1,2,3,4,5]
    plt.figure(figsize=(10,4))
    for part_name, trait_name, param_range, color in zip(part_names[ip], trait_names[ip], param_ranges[ip], colors[ip]):
        print trait_name, (param_range[0], param_range[-1])
        param_arr = np.zeros(len(param_range))
        metric_arr = np.zeros((len(param_range),2))
        for i, trait_val in enumerate(param_range): 
            print i,
            if trait_name == 'R':
                R = trait_val
                modified_plant = plant
            else: 
                R = Rmax
                modified_plant = initialize_modified_plant(sp, soil_type, part_name, trait_name, trait_val)
            gam,eta,k,sw,sst,sCrit=modified_plant.get_derived_params(VPD, s, alpha, n, Ks, sfc, plc=plc)
            ps, assm = simulate_ps_t(n_trajectories, tRun, dt, s0, lam, gam, eta, k, sw, sst, s1, Amax, R)
            if params=='sigma': param_arr[i] = modified_plant.get_sigma(VPD, np.mean(ps, axis=0)[::10])
            elif params=='beta': param_arr[i] = modified_plant.get_Pg12() / modified_plant.stem.P50_stem
            elif params=='alpha': param_arr[i] = modified_plant.canopy.A_canopy/ (modified_plant.soil_root.L_root* modified_plant.soil_root.A_root)  
            elif params=='kappa': param_arr[i] = modified_plant.canopy.Gs_leaf / (0.001*modified_plant.stem.Ksat_stem)
            elif params=='phi': param_arr[i] = R/Amax
            # find performance metrics
            metric_arr[i,0] = get_psCrit(ps, sCrit)
            metric_arr[i,1] = get_relGnet(assm, R)
#             fluxes, _ = modified_plant.get_fluxes(VPD, np.mean(ps, axis=0)[::10])
#             plt.plot(fluxes[:,0])
#             plt.plot(np.mean(ps, axis=0)[::10], label=i)
#             plt.hlines(modified_plant.get_sCrit(plc), 0, tmax)
        plt.subplot(1,2,1); plt.plot(param_arr, metric_arr[:,0], 'o', color=color, label=trait_name); plt.ylabel('Hydraulic risk'); plt.xlabel(params)
        plt.subplot(1,2,2); plt.plot(param_arr, metric_arr[:,1], 'o', color=color, label=trait_name); plt.ylabel('Carbon assimilation'); plt.xlabel(params)
        plt.legend()
        plt.tight_layout()

plot_metrics_param_dependence(juni_plant, VPD, tmax, params='phi')
plot_metrics_param_dependence(pine_plant, VPD, tmax, params='phi')
plt.show()