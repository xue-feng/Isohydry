import matplotlib.pyplot as plt
import cPickle as pickle
import params_constants
import numpy as np
from scipy import stats
from params_soil import soil_dict
from SALib.analyze import sobol
from utility_functions import import_traits_data, get_part
import matplotlib.ticker as mtick

def plot_ind(ax, ind, out, ylims, title=None): 
    ax.plot(ind, out, 'o', color='gray', alpha=0.1, ms=5, mec='None'); # ax.set_title(title)
    spcorr, corr_pval = stats.spearmanr(ind,out); print spcorr, corr_pval
#     if corr_pval<0.05:
#         ax.annotate(r'$r_{s}=$'+str(np.round(spcorr,3))+', '+r'$p=$%.2E' %corr_pval, xy=(0.02,0.9), xycoords='axes fraction', fontsize=8)
    ax.set_ylim(ylims[0], ylims[1])
    
def plot_summary_ind(ax, x,y,xerr,yerr ):
    ax.errorbar(x, y, yerr=yerr, fmt='o', mec='None')
       
def get_nondimen_indices(sp, VPD, tmax):
    with open('../Si_'+sp+'_vpd'+str(int(VPD))+'_tmax'+str(tmax)+'_severeM2M2_outcomes.pickle', 'rb') as handle:
        Y = pickle.load(handle)
    with open('../Si_'+sp+'_vpd'+str(int(VPD))+'_tmax'+str(tmax)+'_severeM2M2_params.pickle', 'rb') as handle:
        params = pickle.load(handle)
    var_dict = {'A_canopy':0,'Gs_leaf':1,'c_leaf':2,'L_stem':3,'A_stem':4,'Ksat_stem':5,'a_stem':6,'P50_stem':7,'L_root':8,'A_root':9,'Rmax':10}
    
    ''' prepare indices '''
    rhoH2O, Mw = params_constants.rhoH2O, params_constants.Mw
    nu, g = params_constants.nu,  params_constants.g
    Sd = params_constants.Sd
    d_root = 0.0005
    gc_max = params[:,var_dict['Gs_leaf']]*params[:,var_dict['A_canopy']]*Mw*Sd/(rhoH2O)
    gx_max = params[:,var_dict['Ksat_stem']]*Sd*params[:,var_dict['A_stem']]/(params[:,var_dict['L_stem']]*rhoH2O)
    gs_max = nu*Ks/(rhoH2O*g)*np.sqrt(params[:,var_dict['A_root']]/(d_root*params[:,var_dict['L_root']]))
    
    beta = params[:, var_dict['P50_stem']]/ (-2.30259/params[:, var_dict['c_leaf']])# convert from c_leaf to P10_leaf
    delta = params[:, var_dict['P50_stem']]/Ps
    kappa = gx_max/gs_max; 
    chi = (gx_max*params[:, var_dict['P50_stem']])/(gc_max*VPD)
    rho = params[:, var_dict['Rmax']]
    epsilon = gs_max*VPD/(alpha*lam)
    tau = tmax*n*var_dict['L_root']/(gc_max*VPD)
    
    HR = Y[:,0]
    Assm = Y[:,1]
    return beta, delta, kappa, chi, rho, epsilon, tau, HR, Assm

def bin_nondimen_groups(data, output1, output2, nbins=12):
    bins = np.linspace(np.min(data), np.max(data), nbins)
    digitized = np.digitize(data, bins)
    binned = np.zeros((nbins,6))
    for i in range(nbins):
        binned[i,0] = np.mean(data[digitized==(i+1)])
        binned[i,1] = np.std(data[digitized==(i+1)])
        
        binned[i,3] = stats.sem(output1[digitized==(i+1)])
        if np.isnan(binned[i,3]): binned[i,2] = np.nan
        else: binned[i,2] = np.mean(output1[digitized==(i+1)])
        
        binned[i,5] = stats.sem(output2[digitized==(i+1)])
        if np.isnan(binned[i,5]): binned[i,4] = np.nan
        else: binned[i,4] = np.mean(output2[digitized==(i+1)])
    return binned

def plot_nondimen_samples(sp='JUNI', VPD=2.0, tmax=180, option='HR'):
    beta, delta, kappa, chi, rho, epsilon, tau, HR, Assm = get_nondimen_indices(sp, VPD, tmax)
    beta_bin = bin_nondimen_groups(beta, HR, Assm)
    delta_bin = bin_nondimen_groups(delta, HR, Assm)
    kappa_bin = bin_nondimen_groups(kappa, HR, Assm)
    chi_bin = bin_nondimen_groups(chi, HR, Assm)
    rho_bin = bin_nondimen_groups(rho, HR, Assm)
    epsilon_bin = bin_nondimen_groups(epsilon, HR, Assm)
    tau_bin = bin_nondimen_groups(tau, HR, Assm)
    
    ''' plotting sample results '''
    if option=='CA': out = Assm; group_ind =(4,5); suptitle = 'C assimilation: %s at %s kPa and %s days'%(sp, VPD, tmax); ylims = (-0.03, 0.12)
    elif option=='HR': out = HR; group_ind =(2,3); suptitle = 'Hydraulic risk: %s at %s kPa and %s days'%(sp, VPD, tmax); ylims=(-0.1,1.0)
    
    fig = plt.figure(figsize=(15,3));# plt.suptitle(suptitle)
    get_inputs = lambda databin: (databin[:,0], databin[:,group_ind[0]], databin[:,1], databin[:,group_ind[1]])
    ax = fig.add_subplot(171)
    plot_ind(ax,beta,out,ylims,'P50/Pg12'); plot_summary_ind(ax,*get_inputs(beta_bin))
    
    ax = fig.add_subplot(172)
    plot_ind(ax,delta,out,ylims,'P50/Psat_soil'); plot_summary_ind(ax,*get_inputs(delta_bin))
    ax.yaxis.set_ticklabels([])
    _, labels = plt.xticks(); plt.setp(labels, rotation=45)
#     ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    
    ax = fig.add_subplot(173)
    plot_ind(ax,kappa,out,ylims,'gmax_stem/gmax_root'); plot_summary_ind(ax,*get_inputs(kappa_bin))
    ax.yaxis.set_ticklabels([])
    _, labels = plt.xticks(); plt.setp(labels, rotation=45)
    
    ax = fig.add_subplot(174)
    plot_ind(ax,chi,out,ylims,'(gmax_stem*P50)/(gmax_canopy*VPD)'); plot_summary_ind(ax,*get_inputs(chi_bin))
    ax.yaxis.set_ticklabels([])
    _, labels = plt.xticks(); plt.setp(labels, rotation=45)
    
    ax = fig.add_subplot(175)
    plot_ind(ax,rho,out,ylims,'Rmax/Amax'); plot_summary_ind(ax,*get_inputs(rho_bin))
    ax.yaxis.set_ticklabels([])
    _, labels = plt.xticks(); plt.setp(labels, rotation=45)
    
    ax = fig.add_subplot(176)
    plot_ind(ax,epsilon,out,ylims,'(gmax_canopy*VPD)/(alpha*lam)'); plot_summary_ind(ax,*get_inputs(epsilon_bin))
    ax.yaxis.set_ticklabels([])
    _, labels = plt.xticks(); plt.setp(labels, rotation=45)
    
    ax = fig.add_subplot(177)
    plot_ind(ax,tau,out,ylims,'tmax*n*Zr/(gmax_canopy*VPD)'); plot_summary_ind(ax,*get_inputs(tau_bin))
    ax.yaxis.set_ticklabels([])
    _, labels = plt.xticks(); plt.setp(labels, rotation=45)
    
#     plt.tight_layout()
    
def define_problem(sp='JUNI'):
    var_vals = [traits[sp][get_part(v)][v] for v in var_names[:-1]]; var_vals.extend([Rmax])
    problem_bounds = [[min(0.5*v,2.0*v), max(0.5*v,2.0*v)] for v in var_vals]
    problem = {'num_vars': n_vars,'names': var_names,'bounds': problem_bounds}
    return problem

def sobol_analysis(sp, VPD, tmax):
    # retrieve sample
    with open('../Si_'+sp+'_vpd'+str(int(VPD))+'_tmax'+str(tmax)+'_severeM2M2_outcomes.pickle', 'rb') as handle:
        Y = pickle.load(handle)
    # perform analysis
    problem = define_problem(sp)
    Si_A = sobol.analyze(problem, Y[:,1], calc_second_order=False, print_to_console=True)
    try: Si_H = sobol.analyze(problem, Y[:,0], calc_second_order=False, print_to_console=True)
    except: pass

    Si_depo = np.zeros((n_vars, 4))
    Si_depo[:,0] = Si_H['ST']
    Si_depo[:,1] = Si_H['ST_conf']
    Si_depo[:,2] = Si_A['ST']
    Si_depo[:,3] = Si_A['ST_conf']
    return Si_depo
 
def verify(sp='JUNI', VPD=2.0, tmax=180):
    Si_depo = sobol_analysis(sp, VPD, tmax)      
    plt.figure(figsize=(6,4.5))
    plt.subplot(211)
    plt.bar(np.arange(len(Si_depo)),Si_depo[:,0], yerr=Si_depo[:,1]); plt.ylim(0,1)
    plt.subplot(212)
    plt.bar(np.arange(len(Si_depo)),Si_depo[:,2], yerr=Si_depo[:,3]); plt.ylim(0,1)
    plt.tight_layout()

def barplot_earlylate(sp='JUNI', VPD=2.0, early=30, late=180):
    Si_early = sobol_analysis(sp, VPD, early)
    Si_late = sobol_analysis(sp, VPD, late)
    
    ind = np.arange(len(Si_early)); width = 0.35
    error_kw = dict(ecolor='darkgray', lw=1, capsize=2, capthick=1)
    
    plt.figure(figsize=(6,5))
    ax = plt.subplot(211)
    ax.bar(ind, Si_early[:,0], width, color='blue', edgecolor='none', yerr=Si_early[:,1], error_kw=error_kw)
    ax.bar(ind+width, Si_late[:,0], width, color='grey', edgecolor='none', yerr=Si_late[:,1], error_kw=error_kw)
    ax.set_xticks(ind + width)
    ax.set_xticklabels([])
    plt.ylim(0,1.2); plt.xlim(-width, np.max(ind)+width*3.0)
    ax = plt.subplot(212)
    ax.bar(ind, Si_early[:,2], width, color='blue', edgecolor='none', yerr=Si_early[:,3], error_kw=error_kw)
    ax.bar(ind+width, Si_late[:,2], width, color='grey', edgecolor='none', yerr=Si_late[:,3], error_kw=error_kw)
    ax.set_xticks(ind + width)
    ax.set_xticklabels(var_names, rotation=45)
    plt.ylim(0,1.2); plt.xlim(-width, np.max(ind)+width*3.0)
    plt.tight_layout()

def barplot_twospecies(VPD=2.0, tmax=30):
    Si_juni = sobol_analysis('JUNI', VPD, tmax)
    Si_pine = sobol_analysis('PINE', VPD, tmax)
    
#     HR = Si_pine[:,0]+Si_juni[:,0]; HRisort = HR.argsort()
#     CA = Si_pine[:,2]+Si_juni[:,2]; CAisort = CA.argsort()
#     iplot = np.argsort(CAisort.argsort() + HRisort.argsort())[::-1][:7]
#     var_names = np.array(['A_canopy','Gs_leaf','c_leaf','L_stem','A_stem','Ksat_stem','a_stem','P50_stem','L_root','A_root','Rmax'])
    iplot = np.array([7,2,0,1,8,9,10])
    
    ind = np.arange(len(iplot)); width = 0.35
    error_kw = dict(ecolor='darkgray', lw=1, capsize=2, capthick=1)
    
    plt.figure(figsize=(5,5))
    ''' for hydraulic risk'''
    ax = plt.subplot(211)
    ax.bar(ind, Si_juni[iplot,0], width, color='brown', edgecolor='white', yerr=Si_juni[iplot,1], error_kw=error_kw)
    ax.bar(ind+width, Si_pine[iplot,0], width, color='blue', alpha=0.5, edgecolor='white', yerr=Si_pine[iplot,1], error_kw=error_kw)
    ax.set_xticks(ind + width)
    ax.set_xticklabels([])
    plt.ylim(0,1.2); plt.xlim(-width, np.max(ind)+width*3.0)
    
    ''' for carbon '''
    ax = plt.subplot(212) 
    ax.bar(ind, Si_juni[iplot,2], width, color='brown', edgecolor='white', yerr=Si_juni[iplot,3], error_kw=error_kw)
    ax.bar(ind+width, Si_pine[iplot,2], width, color='blue', alpha=0.5, edgecolor='white', yerr=Si_pine[iplot,3], error_kw=error_kw)
    ax.set_xticks(ind + width)
    ax.set_xticklabels(var_names[iplot], rotation=45)
    print var_names[iplot]
    plt.ylim(0,1.2); plt.xlim(-width, np.max(ind)+width*3.0)
    plt.tight_layout()
    
## soil conditions ##
soil_type = 'loamy_sand'
smin = soil_dict[soil_type]['sh']
sfc = soil_dict[soil_type]['sfc']
sst = soil_dict[soil_type]['sst']
n = soil_dict[soil_type]['n']
Ks = soil_dict[soil_type]['Ksat']
Ps = soil_dict[soil_type]['Ps']
traits = import_traits_data()

n_trajectories = 500; dt = 0.1
lam=0.05; alpha=0.010;  Amax = 1.0/dt; Rmax = 0.10*Amax;

## prep for sensitivity analysis ##
var_names = np.array(['A_canopy','Gs_leaf','c_leaf','L_stem','A_stem','Ksat_stem','a_stem','P50_stem','L_root','A_root','Rmax'])
n_vars = len(var_names)

VPD=2.0

barplot_twospecies(tmax=30)
barplot_twospecies(tmax=180)
plt.show()
 
# barplot_earlylate('PINE')
# barplot_earlylate('JUNI')
# plt.show()

tmax = 180
plot_nondimen_samples('JUNI', VPD=VPD, tmax=tmax, option='HR')
plot_nondimen_samples('JUNI', VPD=VPD, tmax=tmax, option='CA')
plot_nondimen_samples('PINE', VPD=VPD, tmax=tmax, option='HR')
plot_nondimen_samples('PINE', VPD=VPD, tmax=tmax, option='CA')
plt.show()
