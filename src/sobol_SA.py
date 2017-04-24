from utility_functions import import_traits_data, simulate_ps_t, get_part, initialize_generic_plant
from SALib.sample import saltelli
from SALib.analyze import sobol
from params_soil import soil_dict
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import cPickle as pickle
import time

def get_psCrit(ps, sCrit):
    return len(ps[ps<sCrit])/float(np.shape(ps)[0]*np.shape(ps)[1])

def get_relGnet(assm, R):
    return np.mean(np.sum(assm, axis=1))/((Amax-R)*np.shape(assm)[1]) 

def plot_sensitivity(Si):
    fig = plt.figure(figsize=(7,3))
    ax = fig.add_subplot(111)
    width=0.4
    ax.bar(np.arange(n_vars)-width*0.5,Si['ST'], width, color='gray', yerr=Si['ST_conf']*0.5)
    ax.set_xticks(np.arange(10))
    ax.set_xticklabels(var_names)
    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=30, fontsize=10)
    # plt.ylim(0,1.1)
    plt.xlim(-1,n_vars)
    
def evaluate_model(args):
    param_values, VPD, tmax = args
    Y = np.zeros((len(param_values),2))
    tRun = np.arange(0,tmax+dt,dt)
    for i, p_vals in enumerate(param_values):
        plant_traits = p_vals[:-1]; R = p_vals[-1]
        generic_plant = initialize_generic_plant(var_names, plant_traits, soil_type)
        gam,eta,k,sw,sst,sCrit=generic_plant.get_derived_params(VPD, s, alpha, n, Ks, sfc, plc=plc)
        ps, assm = simulate_ps_t(n_trajectories, tRun, dt, s0, lam, gam, eta, k, sw, sst, s1, Amax, R)
        Y[i,0] = get_psCrit(ps, sCrit)
        Y[i,1] = get_relGnet(assm, R)  
        print i, Y[i]
    return Y

def Sobol(sp, n_runs):
    var_vals = [traits[sp][get_part(v)][v] for v in var_names[:-1]]; var_vals.extend([Rmax])
    problem_bounds = [[min(0.5*v,2.0*v), max(0.5*v,2.0*v)] for v in var_vals]
    problem = {'num_vars': n_vars,'names': var_names,'bounds': problem_bounds}
    # generate samples
    param_values = saltelli.sample(problem, n_runs, calc_second_order=False); print len(param_values)
    
    VPD=2.0; gridsize = 10
    Si_depo = np.zeros((gridsize, n_vars, 4)) # number of increments, parameters, Assm/HF/CIs
    for i, tmax in enumerate(np.linspace(60,240,gridsize)): 
        print i,
        ''' parallel processing '''
        ncores = mp.cpu_count(); print('There are %s cores on this machine '%(str(ncores),))
        pool = mp.Pool()
        param_cores = np.array_split(param_values, ncores, axis=0)
        param_augmented = [(p, VPD, tmax) for p in param_cores]
        Y_pooled = pool.map(evaluate_model, param_augmented)
        Y = np.vstack(Y_pooled) # output shape (len(param_values), 2)
    
        # perform analysis
        Si_A = sobol.analyze(problem, Y[:,1], calc_second_order=False, print_to_console=True)
        try: Si_H = sobol.analyze(problem, Y[:,0], calc_second_order=False, print_to_console=True)
        except: pass
        
        # store 
        Si_depo[i,:,0] = Si_H['ST']
        Si_depo[i,:,1] = Si_A['ST_conf']
        Si_depo[i,:,2] = Si_H['ST']
        Si_depo[i,:,3] = Si_A['ST_conf']
        with open('/Users/xuefeng/Dropbox/Projects/Isohydricity/Sobol/Si_'+sp+'_vpd'+str(int(VPD))+'_tmax'+str(tmax)+'.pickle', 'wb') as handle:
            pickle.dump(Si_depo, handle)
    return Si_depo

# set soil conditions 
soil_type = 'loamy_sand'
smin = soil_dict[soil_type]['sh']
sfc = soil_dict[soil_type]['sfc']
sst = soil_dict[soil_type]['sst']
n = soil_dict[soil_type]['n']
Ks = soil_dict[soil_type]['Ksat']
s = np.hstack((np.linspace(smin,sst,500), np.linspace(sst+0.001,1.0,100)))
traits = import_traits_data()

n_trajectories = 500; dt = 0.1
lam=0.15; alpha=0.010; s1 = sfc; s0 = 0.5; Amax = 1.0/dt; Rmax = 0.10*Amax; plc = 0.5

''' prep for Sobol sensitivity analysis '''
var_names = np.array(['A_canopy','Gs_leaf','c_leaf','L_stem','A_stem','Ksat_stem','a_stem','P50_stem','L_root','A_root','Rmax'])
n_vars = len(var_names)
tmax=240.0; VPD = 2.0

with open('/Users/xuefeng/Dropbox/Projects/Isohydricity/Sobol/Si_JUNI_vpd'+str(int(VPD))+'_tmax'+str(tmax)+'.pickle', 'rb') as handle:
    Si_juni = pickle.load(handle)
with open('/Users/xuefeng/Dropbox/Projects/Isohydricity/Sobol/Si_PINE_vpd'+str(int(VPD))+'_tmax'+str(tmax)+'.pickle', 'rb') as handle:
    Si_pine = pickle.load(handle)
 
def plot_plant_sensitivities(Si, mode):
    ordered = np.argsort(Si[-1,:,mode_dict[mode]])[::-1] # ordered by values at last increment - big to small
    gridsize=10; tRun = np.linspace(60,240,gridsize)
    for i,j in enumerate(ordered):
        if i<5: plt.plot(tRun, Si[:,j,mode_dict[mode]], label=var_names[j])
        else: plt.plot(tRun, Si[:,j,mode_dict[mode]])
    plt.ylabel(mode)
    plt.xlabel('duration (days)')
    plt.legend(loc=2)
      
mode_dict = {'hydraulic':0, 'carbon':1}
mode = 'hydraulic'
plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
plot_plant_sensitivities(Si_juni, mode)
plt.subplot(1,2,2)
plot_plant_sensitivities(Si_pine, mode)
plt.tight_layout()
plt.show()


''' perform Sobol analysis '''
t0=time.time()
Si_depo = Sobol('JUNI',n_runs=100)
plt.figure()
for j,lab in enumerate(var_names):
    plt.plot(Si_depo[:,j,1], label=lab)
plt.legend(loc=2)
plt.figure()
for j,lab in enumerate(var_names):
    plt.plot(Si_depo[:,j,0], label=lab)
plt.legend(loc=2)

Si_depo = Sobol('PINE',n_runs=100)
plt.figure()
for j,lab in enumerate(var_names):
    plt.plot(Si_depo[:,j,1], label=lab)
plt.legend(loc=2)
plt.figure()
for j,lab in enumerate(var_names):
    plt.plot(Si_depo[:,j,0], label=lab)
plt.legend(loc=2)

t1=time.time(); print (t1-t0)/60.0
plt.show()


# # plotting 
# plot_sensitivity(Si_A)
# try: plot_sensitivity(Si_H)
# except: pass 
# plt.show()
