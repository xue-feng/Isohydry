import xlrd
import numpy as np
from plant_soil_model import Canopy, Stem, Soil_root, Whole_plant, Soil

def import_traits_data():
    traits_path = '../hydraulic_traits.xls'
    filepath = traits_path
    book = xlrd.open_workbook(filepath)
    sheet = book.sheet_by_name('parameters')
    keys = np.asarray(filter(None, sheet.col_values(0)), dtype='str')
    species = ['JUNI', 'PINE']; sp_coln = [2,4]
    canopy_keys = ['A_canopy', 'Gs_leaf', 'c_leaf', 'Amax', 'rho']
    stem_keys = ['L_stem', 'A_stem', 'Ksat_stem', 'a_stem', 'plc_form', 'P50_stem']
    root_keys = ['L_root', 'A_root', 'd_root']
    chap_dict = {}
    for sp, spc in zip(species, sp_coln):
        chap_dict[sp] = {}
        for part_keys, part_dict in zip([canopy_keys, stem_keys, root_keys], ['canopy_dict', 'stem_dict', 'root_dict']): 
            chap_dict[sp][part_dict] =  {}
            for key in part_keys:
                j = np.where(keys==key)[0][0]+1         # to account for the column of the species
                chap_dict[sp][part_dict][key] = sheet.col_values(spc)[j] 
    return chap_dict

def initialize_plant(sp, params, soil_type):
    # initializing each species with its name (in strings) and vapor pressure deficit
    canopy_dict = params[sp]['canopy_dict']
    stem_dict = params[sp]['stem_dict']
    root_dict = params[sp]['root_dict']
    plant = Whole_plant(species=sp)
    plant.canopy = Canopy(**canopy_dict)
    plant.stem = Stem(**stem_dict)
    plant.soil_root = Soil_root(soil_type=soil_type, **root_dict)
    plant.soil = Soil(soil_type)
    return plant

def get_part(var):
    if var in ['A_canopy', 'Gs_leaf', 'c_leaf', 'Amax', 'rho']: return 'canopy_dict'
    elif var in ['L_stem','A_stem','Ksat_stem','a_stem','P50_stem']: return 'stem_dict'
    elif var in ['L_root','A_root']: return 'root_dict'
    
def initialize_generic_plant(trait_names, params, soil_type):
    vals = lambda var: params[np.where(var==trait_names)[0][0]]
    sp = 'generic'
    canopy_dict = {'A_canopy':vals('A_canopy'),'Gs_leaf':vals('Gs_leaf'),'c_leaf':vals('c_leaf'), 'Amax':vals('Amax'), 'rho':vals('rho')}
    stem_dict = {'L_stem':vals('L_stem'), 'A_stem':vals('A_stem'), 'Ksat_stem':vals('Ksat_stem'), 'a_stem':vals('a_stem'), 'P50_stem':vals('P50_stem')}
    root_dict = {'L_root':vals('L_root'), 'A_root':vals('A_root'), 'd_root':0.0005}
    
    plant = Whole_plant(species=sp)
    plant.canopy = Canopy(**canopy_dict)
    plant.stem = Stem(**stem_dict)
    plant.soil_root = Soil_root(soil_type=soil_type, **root_dict)
    plant.soil = Soil(soil_type)
    return plant
                
def rho(s, lam, gam, eta, k, sw, sst, s1, Amax, R):
    if s<=sw: rho_s = 0.0; assm_s = 0.0-R
    elif (s>sw)&(s<=sst): rho_s = eta*(s-sw)/(sst-sw); assm_s = Amax*(s-sw)/(sst-sw) - R
    elif (s>sst)&(s<=s1): rho_s = eta; assm_s = Amax-R
    elif (s>s1)&(s<=1.0): rho_s = eta + k*(s-s1)/(1.0-s1); assm_s = Amax-R
    return rho_s, assm_s

def simulate_s_t(depths, tRun, dt, sInit, lam, gam, eta, k, sw, sst, s1, Amax, R):
    ''' simulate for a single trajectory
     to renormalize to per plant basis,  gamma need to be modified!  - matters both for loss and for rain pulse input '''
    s_t = np.zeros(len(tRun))
    assm_t = np.zeros_like(s_t)
    s0 = sInit
    assm_t[0] = 0
    for i in range(len(tRun)): 
        R_normed = depths[i]
        Infil_normed = min(R_normed, 1.0-s0)
        ET_L_normed, ASM = rho(s0, lam, gam, eta, k, sw, sst, s1, Amax, R)
        s_out = max(s0 + Infil_normed - dt*(ET_L_normed), sw)
        # update to next step
        s_t[i] = s_out; s0 = s_out
        assm_t[i] = ASM*dt
    return s_t, assm_t

def simulate_s_t_nonlinearized(depths, tRun, dt, sInit, plant, VPD):
    ''' simulate soil moisture and assimilation trajectories 
    without linearization of E-s and A-s relationships '''
    s_t = np.zeros(len(tRun)); s0 = sInit
    assm_t = np.zeros_like(s_t); assm_t[0] = 0
    R = plant.canopy.R()
    for i in range(len(tRun)): 
        R_normed = depths[i]
        Infil_normed = min(R_normed, 1.0-s0)
        # 1. turn s into psi_s, calculate E from plant, A from gs
        Fluxes, _ = plant.get_fluxes_scalar(VPD,s0); P_leaf=Fluxes[2]
        Ar, Zr, n = plant.soil_root.A_root, plant.soil_root.L_root, plant.soil.n
        Ksat, b = plant.soil.Ksat_soil, plant.soil.b
        L = Ksat*s0**(2*b+3)/(n*Zr*Ar)
        E = Fluxes[0]/(n*Zr*Ar)
        
        ASM = plant.canopy.A(P_leaf)
        s_out = max(s0 + Infil_normed - dt*(E + L), 0)
        # update to next step
        s_t[i] = s_out; s0 = s_out
        assm_t[i] = (ASM-R)*dt
    return s_t, assm_t
    

def simulate_rainfall(n_trajectories, tRun, dt, lam, gam):
    size = len(tRun)*n_trajectories 
    depthExp = -np.log(1.0-np.random.random(size=size))/gam
    freqUnif = np.random.random(size=size)
    
    depth = np.zeros(size)
    # the occurence of rainfall in any independent interval is lam*dt
    depth[freqUnif<np.tile(lam,size)*dt] = depthExp[freqUnif<np.tile(lam,size)*dt] # rain falls according to prob within an increment
    depth_re = np.reshape(depth, (n_trajectories, len(tRun)))
    return depth_re

def simulate_ps_t(n_trajectories, tRun, dt, s0, lam, gam, eta, k, sw, sst, s1, Amax, R):
    ''' simulate for multiple trajectories over time '''
    depth_re = simulate_rainfall(n_trajectories, tRun, dt, lam, gam)
    ps_samples = np.zeros((n_trajectories, len(tRun)))
    assm_samples = np.zeros_like(ps_samples)
    for nsim in range(n_trajectories):
        s_t, assm_t = simulate_s_t(depth_re[nsim],tRun,dt,s0,lam, gam, eta, k, sw, sst, s1, Amax, R)
        ps_samples[nsim] = s_t
        assm_samples[nsim] = assm_t
    return ps_samples, assm_samples

def simulate_ps_t_nonlinearized(n_trajectories, tRun, dt, s0, plant, VPD, lam, alpha):
    Ar, Zr, n = plant.soil_root.A_root, plant.soil_root.L_root, plant.soil.n
    gam = (n*Zr*Ar)/(alpha*Ar)
    depth_re = simulate_rainfall(n_trajectories, tRun, dt, lam, gam)
    ps_samples = np.zeros((n_trajectories, len(tRun)))
    assm_samples = np.zeros_like(ps_samples)
    for nsim in range(n_trajectories):
        s_t, assm_t = simulate_s_t_nonlinearized(depth_re[nsim], tRun, dt, s0, plant, VPD)
        ps_samples[nsim] = s_t
        assm_samples[nsim] = assm_t
    return ps_samples, assm_samples
