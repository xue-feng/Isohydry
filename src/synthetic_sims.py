from utility_functions import import_traits_data, simulate_ps_t, get_part, initialize_generic_plant
from SALib.sample import saltelli
from params_soil import soil_dict
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import cPickle as pickle
import params_constants
from scipy import stats
import time
import sys

# set soil conditions 
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

''' prep for Sobol sensitivity analysis '''
var_names = np.array(['A_canopy','Gs_leaf','c_leaf','L_stem','A_stem','Ksat_stem','a_stem','P50_stem','L_root','A_root','Rmax'])
n_vars = len(var_names)
tmax=240.0; VPD = 2.0
    
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

def generate_samples(sp, n_runs, VPD=2.0, tmax=180):
    var_vals = [traits[sp][get_part(v)][v] for v in var_names[:-1]]; var_vals.extend([Rmax])
    problem_bounds = [[min(0.5*v,2.0*v), max(0.5*v,2.0*v)] for v in var_vals]
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
    
    with open('./Si_'+sp+'_vpd'+str(int(VPD))+'_tmax'+str(tmax)+'_outcomes.pickle', 'wb') as handle:
        pickle.dump(Y, handle)
    with open('./Si_'+sp+'_vpd'+str(int(VPD))+'_tmax'+str(tmax)+'_params.pickle', 'wb') as handle:
        pickle.dump(param_values, handle)
    return Y, param_values

def plot_ind(ind, out, title=None): 
    plt.plot(ind, out, 'o'); plt.title(title)
    slope, intercept, _, p_value, _ = stats.linregress(ind, out)
    predict_y = intercept + slope * ind
    print p_value
    if p_value<0.05: 
        plt.plot(ind, predict_y)

def main():
    sp = 'JUNI'
    
    ''' generate samples '''
    t0=time.time()
    Y, params = generate_samples(sp,n_runs=100, VPD=VPD,tmax=tmax)
    t1=time.time(); print (t1-t0)/60.0
    
#     '''open samples'''
#     sp = 'JUNI'
#     with open('./Si_'+sp+'_vpd'+str(int(VPD))+'_tmax'+str(tmax)+'_outcomes.pickle', 'rb') as handle:
#         Y = pickle.load(handle)
#     with open('./Si_'+sp+'_vpd'+str(int(VPD))+'_tmax'+str(tmax)+'_params.pickle', 'rb') as handle:
#         params = pickle.load(handle)
    var_dict = {'A_canopy':0,'Gs_leaf':1,'c_leaf':2,'L_stem':3,'A_stem':4,'Ksat_stem':5,'a_stem':6,'P50_stem':7,'L_root':8,'A_root':9,'Rmax':10}
     
    beta = (-2.30259/params[:, var_dict['c_leaf']])/params[:, var_dict['P50_stem']] # convert from c_leaf to P10_leaf
    delta = Ps/params[:, var_dict['P50_stem']]
    Sd, rhoH2O, Mw = params_constants.Sd, params_constants.rhoH2O, params_constants.Mw
    nu, rhoH2O, g = params_constants.nu, params_constants.rhoH2O, params_constants.g
    Sd, rhoH2O = params_constants.Sd, params_constants.rhoH2O
    d_root = 0.0005
    gc_max = params[:,var_dict['Gs_leaf']]*params[:,var_dict['A_canopy']]*Mw*Sd/(rhoH2O)
    gx_max = params[:,var_dict['Ksat_stem']]*Sd*params[:,var_dict['A_stem']]/(params[:,var_dict['L_stem']]*rhoH2O)
    gs_max = nu*Ks/(rhoH2O*g)*np.sqrt(params[:,var_dict['A_root']]/(d_root*params[:,var_dict['L_root']]))
    kappa = gs_max/gx_max; 
    f = gc_max*VPD/(gx_max*params[:, var_dict['P50_stem']])
    HR = Y[:,0]
    Assm = Y[:,1]
#     
#     beta_norm = beta/(max(beta)-min(beta))
#     delta_norm = delta/(max(delta)-min(delta))
#     kappa_norm = kappa/(max(kappa)-min(kappa))
#     f_norm = f/(max(f)-min(f))
#     plt.scatter(beta_norm, delta_norm, c=HR, s=50*HR+10, alpha=0.5, lw=0)
#     # plt.figure()
#     # plt.scatter(kappa_norm, f_norm, c=HR, s=50*HR+10, alpha=0.5, lw=0)
#     plt.show()
    
    ''' plotting sample results '''
    fig = plt.figure(figsize=(16,4))
    out = Assm; plt.suptitle('C assimilation')
    # out = HR; plt.suptitle('Hydraulic risk')
    ax = fig.add_subplot(141); plot_ind(beta, out, 'P50/Pg12') 
    ax = fig.add_subplot(142); plot_ind(delta, out, 'Psat_soil/P50')
    ax = fig.add_subplot(143); plot_ind(kappa, out, 'gmax_stem/gmax_root')
    ax = fig.add_subplot(144); plot_ind(f, out, '(gmax_stem*P50)/gmax_canopy*VPD')
    plt.tight_layout()
    plt.savefig('./Si_'+sp+'_vpd'+str(int(VPD))+'_tmax'+str(tmax)+'_Assm.pdf')
    
if __name__ == '__main__':
    main()


