from SALib.sample import saltelli
import numpy as np
from params_soil import soil_dict
import multiprocessing as mp
import cPickle as pickle
import time
from utility_functions import import_traits_data, simulate_ps_t, get_part, initialize_generic_plant
import sys

''' set soil conditions ''' 
soil_type = 'loamy_sand'
smin = soil_dict[soil_type]['sh']
sfc = soil_dict[soil_type]['sfc']
sst = soil_dict[soil_type]['sst']
n = soil_dict[soil_type]['n']
Ks = soil_dict[soil_type]['Ksat']
Ps = soil_dict[soil_type]['Ps']
s = np.hstack((np.linspace(smin,sst,500), np.linspace(sst+0.001,1.0,100)))
traits = import_traits_data()

n_trajectories = 500; dt = 0.1
lam=0.15; alpha=0.010; s1 = sfc; s0 = 0.5; Amax = 1.0/dt; Rmax = 0.10*Amax; plc = 0.5

''' prep for sensitivity analysis '''
var_names = np.array(['A_canopy','Gs_leaf','c_leaf','L_stem','A_stem','Ksat_stem','a_stem','P50_stem','L_root','A_root','Rmax'])
n_vars = len(var_names)
tmax=180; VPD = 2.0; sp='JUNI'; n_runs=10

def get_psCrit(ps, sCrit):
    return len(ps[ps<sCrit])/float(np.shape(ps)[0]*np.shape(ps)[1])

def get_relGnet(assm, R):
    return np.mean(np.sum(assm, axis=1))/((Amax-R)*np.shape(assm)[1]) 

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

def generate_samples(sp, n_runs, VPD=2.0, tmax=180):
    var_vals = [traits[sp][get_part(v)][v] for v in var_names[:-1]]; var_vals.extend([Rmax])
    problem_bounds = [[min(0.5*v,2.0*v), max(0.5*v,2.0*v)] for v in var_vals]
#     problem_bounds = [[min(0.25*v,4.0*v), max(0.25*v,4.0*v)] for v in var_vals]
    problem = {'num_vars': n_vars,'names': var_names,'bounds': problem_bounds}
    # generate samples
    param_values = saltelli.sample(problem, n_runs, calc_second_order=False); print len(param_values)
    
    ''' parallel processing '''
    ncores = mp.cpu_count(); print('There are %s cores on this machine '%(str(ncores),))
    pool = mp.Pool()
    param_cores = np.array_split(param_values, ncores, axis=0)
    param_augmented = [(p, VPD, tmax) for p in param_cores]
    Y_pooled = pool.map(evaluate_model, param_augmented)
    Y = np.vstack(Y_pooled) # output shape (len(param_values), 2)
    return Y, param_values

def main():
    
    ''' generate samples '''
    t0=time.time()
    Y, params = generate_samples(sp,n_runs=n_runs, VPD=VPD,tmax=tmax)
    t1=time.time(); sys.stdout.write('%0.2f minutes \n' %((t1-t0)/60.0) ) 
    
    ''' sample storage '''
    with open('../Si_'+sp+'_vpd'+str(int(VPD))+'_tmax'+str(tmax)+'_outcomes.pickle', 'wb') as handle:
        pickle.dump(Y, handle)
    with open('../Si_'+sp+'_vpd'+str(int(VPD))+'_tmax'+str(tmax)+'_params.pickle', 'wb') as handle:
        pickle.dump(params, handle)
    
if __name__ == '__main__':
    main()
    
    