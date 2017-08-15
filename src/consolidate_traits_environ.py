from utility_functions import import_traits_data, initialize_plant, simulate_ps_t, simulate_ps_t_nonlinearized
from params_soil import soil_dict
import numpy as np
import matplotlib.pyplot as plt
import copy
from scipy import interpolate
import cPickle as pickle
from SALib.sample import saltelli
from SALib.analyze import sobol

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
print 'JUNI Emax', juni_plant.get_Es_params(2.0,s)[-1]
print 'PINE Emax', pine_plant.get_Es_params(2.0,s)[-1]
sigma_juni = juni_plant.get_sigma(2.0,s)
sigma_pine = pine_plant.get_sigma(2.0,s)
print sigma_juni, juni_plant.get_Pg12() - juni_plant.stem.P50_stem
print sigma_pine, pine_plant.get_Pg12() - pine_plant.stem.P50_stem
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

def get_relGnet(assm, Amax, R):
    return np.mean(np.sum(assm, axis=1))/((Amax-R)*np.shape(assm)[1]) 

def simulate_intensity_duration(plant, int_range, dur_range):
    Aref = plant.canopy.gmax_canopy(1.0)
    psCrit_grid = np.zeros((len(int_range), len(dur_range)))
    Gnet_grid = np.zeros_like(psCrit_grid)
    for i, vpd in enumerate(int_range):  # on y-axis of contour
        print i
        for j,tmax in enumerate(dur_range): # on x-axis of contour
            print j,
            # simulation params
            tRun = np.arange(0,tmax+dt,dt)
            # environmental, units in meters
            Avpd = plant.canopy.gmax_canopy(vpd)
            gam,eta,k,sw,sst,sCrit=plant.get_derived_params(vpd, s, alpha, n, Ks, sfc) 
            Amax = plant.canopy.Amax; R = 0.1*Amax
            ps, assm = simulate_ps_t(n_trajectories, tRun, dt, s0, lam, gam, eta, k, sw, sst, s1, Amax, R)
            # probability of being below sCrit over season
            psCrit_grid[i,j] = get_psCrit(ps, sCrit)
            # normalized assm with respect to max net assimilation 
            Gnet_grid[i,j] = get_relGnet(assm, Amax, R)*(Avpd/Aref)
    return psCrit_grid, Gnet_grid

def plot_trajectories(plant, tmax, VPD, newfig=True): 
    if newfig: 
        plt.figure(figsize=(6,8))
    tRun =  np.arange(0,tmax+dt,dt)
    gam, eta, k, sw, sst, sCrit = plant.get_derived_params(VPD, s, alpha, n, Ks, sfc, plc=0.50)
    Amax = plant.canopy.Amax; R = 0.1*Amax
    
#     ps, assm = simulate_ps_t(n_trajectories, tRun, dt, s0, lam, gam, eta, k, sw, sst, s1, Amax, R)
    ps, assm = simulate_ps_t_nonlinearized(n_trajectories, tRun, dt, s0, plant, VPD, lam, alpha)
    
    ps_mean = np.mean(ps, axis=0)
    assm_cumsum = np.cumsum(assm, axis=1)/((Amax-R)*np.shape(assm)[1]*dt)
    
    plt.figure(figsize=(5,5))
    plt.subplot(211); plt.title(str(get_psCrit(ps,sCrit)))
    for i in range(len(ps)):
        plt.plot(tRun, ps[i], lw=0.2, color='darkgray')
        plt.plot(tRun[ps[i]<=sCrit], ps[i][ps[i]<=sCrit], lw=1.5, color='k')
#     plt.plot(tRun, ps_mean, lw=1.5, color='k')
    plt.plot(tRun, np.ones(len(tRun))*sst, color='red', ls=':')
    plt.plot(tRun, np.ones(len(tRun))*sCrit, color='red', ls='--')
    plt.ylabel('Soil moisture')
    plt.subplot(212); plt.title(str(get_relGnet(assm, Amax, R)))
    for i in range(len(ps)):
        plt.plot(tRun, assm_cumsum[i], lw=0.2, color='darkgray')
    plt.plot(tRun, np.mean(assm_cumsum,axis=0), lw=1.5, color='k')
    plt.ylabel('Normalized net assimilation')
    plt.xlabel('Days')
#     plt.ylim(0,0.5)
    plt.tight_layout()
    plt.show()

    fluxes, _ = plant.get_fluxes(VPD, ps_mean)
    plt.subplot(3,1,1); plt.title('Mean soil moisture')
    plt.plot(tRun, ps_mean, lw=1.0, color='gray'); plt.hlines(sCrit,0,tmax, lw=1.0, color='red')
    plt.subplot(3,1,2); plt.title('canopy conductance')
    plt.plot(tRun, plant.canopy.g_canopy(fluxes[:,2], VPD), lw=1.0, color='green')
    plt.subplot(3,1,3); plt.title('xylem conductance')
    plt.plot(tRun, plant.stem.g_stem(fluxes[:,1]), lw=1.0, color='blue')

def plot_intensity_duration(plant, int_range, dur_range):
    TMAX, INT = np.meshgrid(dur_range, int_range)
    _, Gnet_grid = simulate_intensity_duration(plant, int_range, dur_range)
     
    plt.figure(figsize=(6,5.5))
    CS = plt.contour(TMAX, INT, Gnet_grid,15)
    plt.clabel(CS, fontsize=10, inline=1)
    plt.xlabel('Duration (days)'); plt.ylabel('Intensity (kPa)')
    plt.tight_layout()
    
def plot_intensity_duration_sensitivity(plant, int_range, dur_range, xf=50, interp_grid=100):
    _, Gnet_grid = simulate_intensity_duration(plant, int_range, dur_range)
    
    x, y = np.meshgrid(dur_range, int_range*xf) # in dPa 
    tck = interpolate.bisplrep(x,y,Gnet_grid,s=0)
    int_range2 = np.linspace(np.min(int_range)*xf,np.max(int_range)*xf,interp_grid) 
    dur_range2 = np.linspace(np.min(dur_range),np.max(dur_range),interp_grid)
    x2, y2 = np.meshgrid(dur_range2, int_range2)
    z2 = interpolate.bisplev(x2[0,:], y2[:,0], tck)
    dzdx = interpolate.bisplev(x2[0,:], y2[:,0], tck, dx=1)
    dzdy = interpolate.bisplev(x2[0,:], y2[:,0], tck, dy=1)
    
    plt.figure(figsize=(6,5.5))
    z2ma = np.ma.masked_where(dzdx>dzdy, z2)
    plt.pcolor(z2ma)
    
def initialize_modified_plant(sp, soil_type, part_name, trait_name, trait_val):
    traits_modified = copy.deepcopy(traits)
    traits_modified[sp][part_name][trait_name] = trait_val
    plant_modified = initialize_plant(sp, traits_modified, soil_type)
    return plant_modified

def get_iso_aniso_performance(sp, VPD, tmax_arr, iso_xf, aniso_xf, plc=0.8):
    plant = initialize_plant(sp, traits, soil_type)
    
    # defining what it means to be iso or aniso
    Pg12_ref = plant.get_Pg12(); Pg12_iso = iso_xf*Pg12_ref; Pg12_aniso = aniso_xf*Pg12_ref
    plant_iso=initialize_modified_plant(sp,soil_type,'canopy_dict','c_leaf',np.log(0.12)/(Pg12_iso))
    plant_aniso = initialize_modified_plant(sp,soil_type,'canopy_dict','c_leaf',np.log(0.12)/(Pg12_aniso)) 
    
    Aref_iso = plant_iso.canopy.gmax_canopy(0.5); Avpd_iso = plant_iso.canopy.gmax_canopy(VPD)
    Aref_aniso = plant_aniso.canopy.gmax_canopy(0.5); Avpd_aniso = plant_aniso.canopy.gmax_canopy(VPD)
     
    metrics_iso = np.zeros((len(tmax_arr),2))
    metrics_aniso = np.zeros_like(metrics_iso)
    print 'start iso_aniso_performance, loop length:', len(tmax_arr)
    for i, tmax in enumerate(tmax_arr): 
        print i,
        tRun = np.arange(0,tmax+dt,dt)
        ## isohydric performance
        gam,eta,k,sw,sst,sCrit=plant_iso.get_derived_params(VPD, s, alpha, n, Ks, sfc, plc)
        Amax = plant.canopy.Amax; R = 0.1*Amax
        ps, assm = simulate_ps_t(n_trajectories, tRun, dt, s0, lam, gam, eta, k, sw, sst, s1, Amax, R)
        metrics_iso[i,0] = get_psCrit(ps, sCrit)
        metrics_iso[i,1] = get_relGnet(assm, Amax, R) *(Avpd_iso/Aref_iso)
#         for k in range(len(ps[:,0])):
#             plt.plot(tRun, ps[k], lw=0.5, color='gray')
#         plt.hlines(sw, 0, tmax); plt.hlines(sst, 0, tmax); plt.hlines(sCrit, 0, tmax, 'red')
#         plt.show()
        ## anisohydric performance
        gam,eta,k,sw,sst,sCrit=plant_aniso.get_derived_params(VPD, s, alpha, n, Ks, sfc, plc)
        Amax = plant.canopy.Amax; R = 0.1*Amax
        ps, assm = simulate_ps_t(n_trajectories, tRun, dt, s0, lam, gam, eta, k, sw, sst, s1, Amax, R)
        metrics_aniso[i,0] = get_psCrit(ps, sCrit)
        metrics_aniso[i,1] = get_relGnet(assm, Amax, R) *(Avpd_aniso/Aref_aniso)
#         for k in range(len(ps[:,0])):
#             plt.plot(tRun, ps[k], lw=0.5, color='gray')
#         plt.hlines(sw, 0, tmax); plt.hlines(sst, 0, tmax); plt.hlines(sCrit, 0, tmax, 'red')
#         plt.show()
    return metrics_iso, metrics_aniso
    
def plot_iso_aniso_performance(sp, VPD, tmax_arr, iso_xf, aniso_xf, plc=0.8):
    ''' plot the hydaulic risk and carbon gain of isohydric and anisohydric modified_plants as function of tmax 
        isohydry and anisohydry defined in terms of shifts in Pg12
    '''
    metrics_iso, metrics_aniso = plot_iso_aniso_performance(sp, VPD, tmax_arr, iso_xf, aniso_xf, plc)
        
    plt.figure(figsize=(4,6))
    plt.subplot(2,1,1)
    plt.plot(tmax_arr, metrics_iso[:,0], label='iso'); plt.plot(tmax_arr, metrics_aniso[:,0], label='aniso')
    plt.ylim(-0.1,1); plt.xlim(0,np.max(tmax_arr)); plt.xlabel('Duration (days)'); plt.ylabel('"Hydraulic risk"')
    plt.legend()
    plt.title(sp+', VPD: '+str(VPD)+'kPa')
    plt.subplot(2,1,2)
    plt.plot(tmax_arr, metrics_iso[:,1]); plt.plot(tmax_arr, metrics_aniso[:,1])
    plt.xlim(0,np.max(tmax_arr)); plt.xlabel('Duration (days)'); plt.ylabel('"Carbon gain"')
    plt.tight_layout()

    
# visualize pine and juniper performance under intensity and duration 
# performance gauged by hydraulic risk and reduction in Assimilation

n_trajectories = 50; dt = 0.1 #dt = 0.1
# lam=0.15; alpha=0.010; s1 = sfc; s0 = 0.5; Amax = 1.0/dt; R = 0.10*Amax
lam=0.05; alpha=0.007; s1 = sfc; s0 = sst; # Amax = 1.0/dt; R = 0.10*Amax
gridsize=10; int_range=np.linspace(0.5,4.0,gridsize); dur_range=np.linspace(30,180,gridsize)

plot_trajectories(juni_plant, 180, 2.0, newfig=True)
# plot_trajectories(pine_plant, 180, 2.0, newfig=True)
plt.show()

# plot_sigma_VPD(np.linspace(0.5,4.0,10))
# plt.show()

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

# gridsize=10; int_range=np.linspace(0,3.0,gridsize); dur_range=np.linspace(0,180,gridsize)
# plot_intensity_duration(pine_plant, int_range, dur_range)
# plot_intensity_duration(juni_plant, int_range, dur_range)
# plt.show()

''' showing regions of more sensitivity to duration or intensity '''
# gridsize=50; int_range=np.linspace(0.5,3.0,gridsize); dur_range=np.linspace(30,180,gridsize)
# plot_intensity_duration_sensitivity(juni_plant, int_range, dur_range, xf=50, interp_grid=50)
# plot_intensity_duration_sensitivity(pine_plant, int_range, dur_range, xf=50, interp_grid=50)
# plt.show()

# def plot_crits(plant, gridsize):
#     int_range=np.linspace(0.5,3.0,gridsize); dur_range=np.linspace(30,180,gridsize)
#     _, Gnet_grid = simulate_intensity_duration(plant, int_range, dur_range)
# 
#     dzdx = np.diff(Gnet_grid,axis=0)[:,1:]
#     dzdy = np.diff(Gnet_grid,axis=1)[1:,:]
#     
#     xc, yc = np.where(dzdx>dzdy)
#     pos = np.where(np.diff(xc))[0]+1
#     xgroups = np.split(xc,pos)
#     ygroups = np.split(yc,pos)
#     
#     ycrit = np.zeros(len(Gnet_grid[0])-1)
#     for xg, yg in zip(xgroups, ygroups):
#         ycrit[xg[0]] = np.max(yg)
#     plt.figure(figsize=(6,5.5))
#     plt.plot(dur_range[1:], np.take(int_range[1:], np.array(ycrit,dtype=int)), 'o')
# 
#     
#     with open('/Users/xuefeng/Desktop/'+plant.species+'_crits.pickle', 'wb') as output_file:
#         pickle.dump([dur_range[1:], np.take(int_range[1:], np.array(ycrit,dtype=int))], output_file)
# 
# plot_crits(juni_plant,50)
# plot_crits(pine_plant,50)
# plt.show()


# dur_range = np.linspace(0,240,10)
# plot_iso_aniso_performance(sp='JUNI', VPD=0.5, tmax_arr=dur_range, iso_xf=0.5, aniso_xf=2.0, plc=0.5)
# plt.show()

# iso, aniso = get_iso_aniso_performance(sp='PINE', VPD=0.5, tmax_arr=dur_range, iso_xf=0.5, aniso_xf=2.0, plc=0.5)
# plt.figure(figsize=(4,6))
# plt.subplot(2,1,1)
# plt.plot(dur_range, iso[:,0], label='iso'); plt.plot(dur_range, aniso[:,0], label='aniso')
# plt.ylim(-0.1,1); plt.xlim(0,np.max(dur_range)); plt.xlabel('Duration (days)'); plt.ylabel('"Hydraulic risk"')
# plt.legend()
# plt.subplot(2,1,2)
# plt.plot(dur_range, iso[:,1]); plt.plot(dur_range, aniso[:,1])
# plt.xlim(0,np.max(dur_range)); plt.xlabel('Duration (days)'); plt.ylabel('"Carbon gain"')
# plt.tight_layout()
# plt.show()

plant = juni_plant
with open('/Users/xuefeng/Desktop/'+plant.species+'_aniso_limit.pickle', 'rb') as output_file:
    tmax_crit_juni, VPD_range = pickle.load(output_file)
plant = pine_plant
with open('/Users/xuefeng/Desktop/'+plant.species+'_aniso_limit.pickle', 'rb') as output_file:
    tmax_crit_pine, VPD_range = pickle.load(output_file)
print 'pause'

plt.figure(figsize=(6,5.5))
plt.plot(tmax_crit_pine, VPD_range, 'o')
plt.plot(tmax_crit_juni, VPD_range, 'o')
plt.xlabel('Critical tmax (days)')
plt.ylabel('VPD')
plt.xlim(10,180)
plt.show()

# def plot_aniso_limits(plant, gsize):
#     dur_range = np.linspace(0,180,gsize)
#     VPD_range = np.linspace(0.1,3.0,gsize)
#     tmax_crit = np.zeros(gsize)
#     for i, vpd in enumerate(VPD_range):
#         print i, vpd
#         iso, aniso = get_iso_aniso_performance(plant.species, VPD=vpd, tmax_arr=dur_range, iso_xf=0.5, aniso_xf=2.0, plc=0.5)
#     #     plt.plot(aniso[:,0] - iso[:,0]); plt.plot(aniso[:,1] - iso[:,1]); plt.show()
#         imorerisk = np.where(aniso[:,0]-iso[:,0]>aniso[:,1] - iso[:,1])[0]
#         if len(imorerisk)==0: 
#             tmax_crit[i] = dur_range[0]
#         else:
#             tmax_crit[i] = dur_range[imorerisk[0]]
#     plt.plot(tmax_crit, VPD_range, 'o')
#     with open('/Users/xuefeng/Desktop/'+plant.species+'_aniso_limit.pickle', 'wb') as output_file:
#         pickle.dump([tmax_crit, VPD_range], output_file)
# 
# plt.figure(figsize=(6,5.5))
# plot_aniso_limits(juni_plant, 20)
# # plot_aniso_limits(pine_plant, 20)
# plt.show()

''' 
check sigma under given traits for JUNI (invalid, but why?, 0.50-0.52 range) and PINE (valid, 0.49-0.55 range)
derive sigma and Pg12-P50 metrics for trait changes -- INTERNAL CONSISTENCY
think about performance implications (different drought regimes), for more species (MAHERALI?)
KEY: INCORPORATE SOIL WATER LEAKAGE EVEN AFTER STOMATAL CLOSURE!!!
GET A FEW MORE SPECIES FOR TRAITS '''

