import numpy as np
import matplotlib.pyplot as plt
from utility_functions import simulate_ps_t

def plot_environment_contours(sInit,lam,gam,eta,k,sw,sst,s1,Amax,R,int_arr,int_option='eta'):
    tmax_arr = np.linspace(30,180,10)
    psCrit_grid = np.zeros((len(int_arr), len(tmax_arr)))
    Gnet_grid = np.zeros_like(psCrit_grid)

    for i, intensity in enumerate(int_arr): # on y-axis of contour
        print i,
        if int_option=='lam': sim_params = (intensity,gam,eta,k,sw,sst,s1,Amax,R); DI = (gam*eta)/intensity
        if int_option=='eta': sim_params = (lam,gam,intensity,k,sw,sst,s1,Amax,R); DI = (gam*intensity)/lam
        print DI
        for j,tmax in enumerate(tmax_arr): # on x-axis of contour
            print j, 
            tRun = np.arange(0,tmax+dt,dt)
            ps_samples, assm_samples = simulate_ps_t(n_trajectories,tRun,dt,sInit,*sim_params)
            # probability of being below sCrit over season
            psCrit_grid[i,j] = len(ps_samples[ps_samples<sCrit])/float(n_trajectories*len(tRun))
            # normalized assm with respect to max net assimilation 
            Gnet_grid[i,j] = np.mean(np.sum(assm_samples, axis=1))/((Amax-R)*tmax)
            
    TMAX, INT = np.meshgrid(tmax_arr, int_arr)
    try:
        plt.figure(figsize=(4.8,4.5))
        CS = plt.contour(TMAX, INT, psCrit_grid)
        plt.clabel(CS, fontsize=10, inline=1)
        if int_option=='lam': plt.gca().invert_yaxis()
    except ValueError:
        print "Zeros in grid." 
        pass
    
    try:
        plt.figure(figsize=(4.8,4.5))
        CS = plt.contour(TMAX, INT, Gnet_grid)
        plt.clabel(CS, fontsize=10, inline=1)
        if int_option=='lam': plt.gca().invert_yaxis()
    except ValueError:
        print "Zeros in grid." 
        pass
    
    try:
        plt.figure(figsize=(4.8,4.5))
        CS = plt.contour(TMAX, INT, Gnet_grid-psCrit_grid)
        plt.clabel(CS, fontsize=10, inline=1)
        if int_option=='lam': plt.gca().invert_yaxis()
    except ValueError:
        print "Zeros in grid." 
        pass
    
def plot_trait_dependence(n_trajectories,tRun,sInit,lam,gam,eta,k,s1,Amax,R):
    # use the range for juniper in sandy loam
    sst_arr = np.linspace(0.25,0.45,5)
    sw_arr = np.linspace(0.15,0.22,5)
    psCrit_grid = np.zeros((len(sst_arr), len(sw_arr)))
    Gnet_grid = np.zeros_like(psCrit_grid)
    tmax = np.max(tRun)
    for i, _ in enumerate(sst_arr):
        print i 
        for j, _ in enumerate(sw_arr):
            print j, 
            sw, sst = sw_arr[j],sst_arr[i]
            sCrit = sw+0.1*(sst-sw)
            sim_params = (lam,gam,eta,k,sw, sst,s1,Amax,R)
            ps_samples, assm_samples = simulate_ps_t(n_trajectories,tRun,sInit,*sim_params)
            psCrit_grid[i,j] = len(ps_samples[ps_samples<sCrit])/float(n_trajectories*len(tRun))
            Gnet_grid[i,j] = np.mean(np.sum(assm_samples, axis=1))/((Amax-R)*tmax)
    SW, SST = np.meshgrid(sw_arr, sst_arr)
    try:
        plt.figure(figsize=(4.8,4.5))
        CS = plt.contour(SW, SST, Gnet_grid-psCrit_grid)
        plt.clabel(CS, fontsize=10, inline=1)
        plt.xlabel('sw'); plt.ylabel('sst')
    except ValueError:
        print "Zeros in grid." 
        pass
    
    try:
        plt.figure(figsize=(4.8,4.5))
        CS = plt.contour(SW, SST, psCrit_grid)
        plt.clabel(CS, fontsize=10, inline=1)
        plt.xlabel('sw'); plt.ylabel('sst')
        print psCrit_grid
    except ValueError:
        print "Zeros in grid." 
        pass
    
    try:
        plt.figure(figsize=(4.8,4.5))
        CS = plt.contour(SW, SST, Gnet_grid)
        plt.clabel(CS, fontsize=10, inline=1)
        plt.xlabel('sw'); plt.ylabel('sst')
        print Gnet_grid
    except ValueError:
        print "Zeros in grid." 
        pass
    
# plt.figure()
# for n in range(10):
#     plt.plot(tRun, ps_samples[n],lw=0.5,color='lightgrey')
#     plt.plot(tRun, assm_samples[n],lw=0.5,color='pink')
# plt.plot(tRun, ps_samples[-1],lw=1,color='grey')
# # plt.plot(tRun, assm_samples[-1],lw=0.5,color='red')
# plt.plot(tRun, np.mean(ps_samples,axis=0),lw=2,color='grey')
# plt.plot(tRun, np.mean(assm_samples,axis=0),lw=2,color='red')
# plt.hlines(sst, 0, tmax, linestyle='--', lw=1, color='red')
# plt.hlines(sCrit, 0, tmax, linestyle='--', lw=1, color='grey')
# plt.title('Short and severe')
# plt.xlabel('days')
# plt.show()

dt = 0.1; tmax=50
tRun = np.arange(0,tmax+dt,dt)
n_trajectories = 1000

sp='JUNI'
if sp=='PINE': Zr = 200; ETmax = 1.05
if sp=='JUNI': Zr = 300; ETmax = 0.45
lam=0.10; alpha=0.7
n=0.45; Ks=20; s1=0.7 # in cm or cm/day - loamy soil? 
gam = n*Zr/alpha; eta = ETmax/(n*Zr); k = Ks/(n*Zr); print 'reference eta (subject to change):', eta
sw=0.08; sst = 0.35; sCrit = sw+0.1*(sst-sw)
sInit=0.45
Amax = 10; R = 0.1*Amax
''' 
1. need to change Amax based on some reference during non-stressed (e.g., low VPD) conditions
2. need to also calculate sCrit based on changes in P50, or Pg12, etc.
3. check what PINE and JUNI contour would look like wrt to intensity and duration on the axes. 
'''
# lam_arr = np.linspace(0.05,0.5,10)
# plot_environment_contours(sInit,lam,gam,eta,k,sw,sst,s1,Amax,R,lam_arr,'lam')
eta_arr = np.linspace(0.01,0.2,10)
plot_environment_contours(sInit,lam,gam,eta,k,sw,sst,s1,Amax,R,eta_arr,'eta')
plt.show()

## tmax = 60,120,180; eta = 0.004, 0.002 # actually transpiration is VPD sensitive - stomata close
# factor: 1.2, 1.3, 1.5, 2.0, 3.0
# tmax=60; eta = 0.01167/3.0; print 'eta realized:', eta
# tmax=60; eta = 0.01167; print 'eta realized:', eta
tmax=90; eta = 0.0033; print 'eta realized:', eta

# tRun = np.arange(0,tmax+dt,dt)
# plot_trait_dependence(n_trajectories,tRun,sInit,lam,gam,eta,k,s1,Amax,R)
plt.show()

'''
1. carbon assimilation shows an interesting nonlinear dynamic under stressed vs. nonstressed conditions 
this is then partially determined by rooting depth! 
2. shape of the contours determined also by the "type" of drought -- whether due to decrease in rainfall or increase in VPD, need to use Dryness Index as a consolidating metric 
3. use carbon cost of refilling to consolidate crossing with carbon? 
4. plot as function of sst and sw (how do these follow from P50, P0leaf, etc.)? 

3a. if roots are deep enough, then with no cost to refilling, it only makes sense to be more risky
'''

# for n in range(n_trajectories):
#     plt.plot(tRun, ps_samples[n],lw=0.5,color='lightgrey')
#     plt.plot(tRun, assm_samples[n],lw=0.5,color='pink')
#     plt.hlines(sst, 0, tmax, lw=0.5, color='red')
# plt.show()