from utility_functions import import_traits_data, initialize_plant, simulate_ps_t, simulate_ps_t_nonlinearized
from params_soil import soil_dict
import numpy as np
import matplotlib.pyplot as plt
import copy

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

sigma_juni = juni_plant.get_sigma(2.0,s)
sigma_pine = pine_plant.get_sigma(2.0,s)
print 'JUNI sigma, pg12-p50', sigma_juni, juni_plant.get_Pg12() - juni_plant.stem.P50_stem
print 'PINE sigma, pg12-p50', sigma_pine, pine_plant.get_Pg12() - pine_plant.stem.P50_stem
print traits['JUNI']['canopy_dict']
print traits['JUNI']['stem_dict']
print traits['JUNI']['root_dict']

def plot_sigma_VPD(VPD_arr = np.linspace(0.5,4.0,20)):
    sigma_arr = np.zeros((len(VPD_arr),2))
    plt.figure(figsize=(5,4.5))
    for i, vpd in enumerate(VPD_arr):
        sigma_arr[i,0] = juni_plant.get_sigma(vpd, s)
        sigma_arr[i,1] = pine_plant.get_sigma(vpd, s)
    plt.figure(figsize=(5,4.5))
    plt.plot(VPD_arr, sigma_arr[:,0], label='Juniper')
    plt.plot(VPD_arr, sigma_arr[:,1], label='Pine')
    plt.legend()
    
def plot_ETmax_VPD(VPD_arr = np.linspace(0.5,4.0,20)):
    Emax_arr = np.zeros((len(VPD_arr),2))
    for i, vpd in enumerate(VPD_arr):
        Emax_arr[i,0] = juni_plant.get_Es_params(vpd, s)[-1]
        Emax_arr[i,1] = pine_plant.get_Es_params(vpd, s)[-1]
    plt.figure(figsize=(6,4.5))
    plt.plot(VPD_arr, Emax_arr[:,0], label='Juniper')
    plt.plot(VPD_arr, Emax_arr[:,1], label='Pine')
    plt.legend(); plt.tight_layout()

def plot_flux_VPD_lines(plant, P_soil_arr, VPD_arr, g_size=10):
    E_grid = np.zeros((len(VPD_arr), len(P_soil_arr)))
    for i, vpd in enumerate(VPD_arr):
        for j, psoil in enumerate(P_soil_arr):
            E_grid[i,j], _ ,_= plant.flux_solver(psoil, vpd)
    for k, psoil in enumerate(P_soil_arr):
        plt.plot(VPD_arr, E_grid[:,k])
    plt.xlabel('VPD (kPa)')
    plt.ylabel('Transpiration')
    plt.tight_layout()
    
def get_psCrit(ps, sCrit):
    return len(ps[ps<sCrit])/float(np.shape(ps)[0]*np.shape(ps)[1])

def get_relGnet(assm, Amax, R, tmax):
    return np.mean(np.sum(assm, axis=1))/((Amax-R)*tmax) 

def plot_trajectories(plant, tmax, VPD, newfig=True, nonlinear=True): 
    if newfig: 
        plt.figure(figsize=(6,8))
    tRun =  np.arange(0,tmax+dt,dt)
    gam, eta, k, sw, sst, sCrit = plant.get_derived_params(VPD, s, alpha, n, Ks, sfc, plc=0.50)
    Amax = plant.canopy.Amax; R = plant.canopy.R
    
    if nonlinear: 
        ps, assm = simulate_ps_t_nonlinearized(n_trajectories, tRun, dt, s0, plant, VPD, lam, alpha)
    else: 
        ps, assm = simulate_ps_t(n_trajectories, tRun, dt, s0, lam, gam, eta, k, sw, sst, s1, Amax, R)

    ps_mean = np.mean(ps, axis=0)
    assm_cumsum = np.cumsum(assm, axis=1)/((Amax-R)*np.shape(assm)[1]*dt)
    
    plt.figure(figsize=(5,5))
    plt.subplot(211); plt.title(str(get_psCrit(ps,sCrit)))
    for i in range(len(ps)):
        plt.plot(tRun, ps[i], lw=0.2, color='darkgray')
        plt.plot(tRun[ps[i]<=sCrit], ps[i][ps[i]<=sCrit], lw=1.5, color='k')

    plt.plot(tRun, np.ones(len(tRun))*sst, color='red', ls=':')
    plt.plot(tRun, np.ones(len(tRun))*sCrit, color='red', ls='--')
    plt.ylabel('Soil moisture')
    plt.subplot(212); plt.title(str(get_relGnet(assm, Amax, R, tmax)))
    for i in range(len(ps)):
        plt.plot(tRun, assm_cumsum[i], lw=0.2, color='darkgray')
    plt.plot(tRun, np.mean(assm_cumsum,axis=0), lw=1.5, color='k')
    plt.ylabel('Normalized net assimilation')
    plt.xlabel('Days')
    plt.tight_layout()
 
    plt.figure(figsize=(6,8))
    fluxes, _ = plant.get_fluxes(VPD, ps_mean)
    plt.subplot(3,1,1); plt.title('Mean soil moisture')
    plt.plot(tRun, ps_mean, lw=1.0, color='gray'); plt.hlines(sCrit,0,tmax, lw=1.0, color='red')
    plt.subplot(3,1,2); plt.title('canopy conductance, m3/d')
    plt.plot(tRun, plant.canopy.g_canopy(fluxes[:,2]), lw=1.0, color='green')
    plt.subplot(3,1,3); plt.title('xylem conductance, m3/d')
    plt.plot(tRun, plant.stem.g_stem(fluxes[:,1]), lw=1.0, color='blue')

    
def initialize_modified_plant(sp, soil_type, part_name, trait_name, trait_val):
    traits_modified = copy.deepcopy(traits)
    traits_modified[sp][part_name][trait_name] = trait_val
    plant_modified = initialize_plant(sp, traits_modified, soil_type)
    return plant_modified

n_trajectories = 50
lam=0.05; alpha=0.007; s1 = sfc; s0 = sst; 
gridsize=10; dt = 0.1 

plot_trajectories(juni_plant, 180, 2.0, newfig=True)
# plot_trajectories(pine_plant, 180, 2.0, newfig=True)
plt.show()

plot_sigma_VPD(np.linspace(0.5,4.0,10))
plt.show()

P_soil_arr=np.linspace(-10,-0.01,1000)
E_grid = np.zeros((len(P_soil_arr),2))
for i, psoil in enumerate(P_soil_arr):
    E_grid[i,0], _ ,_ = pine_plant.flux_solver(psoil, VPD=0.85)
    E_grid[i,1], _ ,_ = juni_plant.flux_solver(psoil, VPD=2.5)
plt.figure(figsize=(5,4.5))
plt.plot(-P_soil_arr, E_grid[:,0])
plt.plot(-P_soil_arr, E_grid[:,1])
plt.xlabel('P_soil (-MPa)')
plt.ylabel('Transpiration')
plt.xscale('log')

plt.figure(figsize=(5,4.5))
plot_flux_VPD_lines(juni_plant, P_soil_arr=np.linspace(-3.0,-0.5,2), VPD_arr=np.linspace(0.5,4.0,50))
plot_flux_VPD_lines(pine_plant, P_soil_arr=np.linspace(-3.0,-0.5,2), VPD_arr=np.linspace(0.5,4.0,50))
plt.show()
