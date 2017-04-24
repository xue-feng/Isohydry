from utility_functions import import_traits_data, initialize_plant
from params_soil import soil_dict
import numpy as np
import matplotlib.pyplot as plt

def P(beta,d,f,k): 
    return 9*s*f - 10*k*beta*s**(2*b)

def Qp(beta,d,f,k):
    return 10*s*f*beta + k*(s**b)*d*(9*f - 10*beta) 

def Rp(beta,d,f,k):
    return -5*(s**b)*beta*k*(d + 2*(s**b)) + f*(9*k*s**(2*b) + s*(5*beta+9))

def Hpp(beta,d,f,k):
    return -9*d*s**(-b+1) + s*(5*beta+9) + 9*k*s**(2*b) + 5*beta*k*(s**b)*(d-2*(s**b))/f

def eps1(beta,d,f,k):
    return -k*s**(2*b-1)/P(beta,d,f,k)*(Hpp(beta,d,f,k) + 1/f*np.sqrt(-2*P(beta,d,f,k)*Qp(beta,d,f,k)+Rp(beta,d,f,k)**2))

def eps2(beta,d,f,k):
    return -k*s**(2*b-1)/P(beta,d,f,k)*(Hpp(beta,d,f,k) - 1/f*np.sqrt(-2*P(beta,d,f,k)*Qp(beta,d,f,k)+Rp(beta,d,f,k)**2))

def tau1(beta,d,f,k):
    return 1.0/P(beta,d,f,k)*(Hpp(beta,d,f,k) + np.sqrt(-2*P(beta,d,f,k)*Qp(beta,d,f,k)+Rp(beta,d,f,k)**2)) 
    
def tau2(beta,d,f,k):
    return 1.0/P(beta,d,f,k)*(Hpp(beta,d,f,k) - np.sqrt(-2*P(beta,d,f,k)*Qp(beta,d,f,k)+Rp(beta,d,f,k)**2))

def get_traits(plant):
    gX = plant.stem.gmax_stem(); psi_T50 = plant.stem.P50_stem
    gS = plant.soil_root.gmax_root(); psi_ssat = Ps
    gC = plant.canopy.gmax_canopy(VPD); psi_L90 = plant.get_Pg12()
    return gX, psi_T50, gS, psi_ssat, gC, psi_L90

def get_params(gX, psi_T50, gS, psi_ssat, gC, psi_L90):
    beta = psi_L90/psi_T50
    d = psi_ssat/psi_T50
    f = gC*vpd/(gX*psi_T50)
    k = gS/gX
    return beta, d, f, k

# set soil conditions 
soil_type = 'loamy_sand'
smin = soil_dict[soil_type]['sh']
sfc = soil_dict[soil_type]['sfc']
sst = soil_dict[soil_type]['sst']
n = soil_dict[soil_type]['n']
Ks = soil_dict[soil_type]['Ksat']
b = soil_dict[soil_type]['b'] 
Ps = soil_dict[soil_type]['Ps'] 
s = 0.2; psi_s = Ps*s**(-b)
VPD=2.0; vpd = VPD/101.325 # VPD in nondimensional form




import cPickle as pickle
vpd_dry=2.0
with open('/Users/xuefeng/Dropbox/Projects/Crossing/Datasets/California_desert/Es_params_VPD'+str(vpd_dry)+'.pickle', 'rb') as handle:
    Es_dict = pickle.load(handle)
sst_arr = np.zeros((len(Es_dict.keys()),2))
for i, sp in enumerate(Es_dict.keys()):
    sst_arr[i,0] = Es_dict[sp][1]
    sst_arr[i,1] = Ps*sst_arr[i,0]**(-b)
    
print 'pause'
# sw_sp, sst_sp, Emax_sp, s_P90, s_Pg12, Zr_sp = Es_dict[sp]


# set species traits
traits = import_traits_data()
sp = 'JUNI'
plant = initialize_plant(sp, traits, soil_type)
gX = plant.stem.gmax_stem(); psi_T50 = plant.stem.P50_stem
gS = plant.soil_root.gmax_root(); psi_ssat = Ps
gC = plant.canopy.gmax_canopy(VPD); psi_L90 = plant.get_Pg12()

gX, psi_T50, gS, psi_ssat, gC, psi_L90 = get_traits(plant)
beta, d, f, k = get_params(gX, psi_T50, gS, psi_ssat, gC, psi_L90)

shift_range = np.linspace(0.5, 2.0, 10)
psi_L90_arr = shift_range*psi_L90
psi_T50_arr = shift_range*psi_T50
beta_arr = np.zeros(len(shift_range))
d_arr = np.zeros_like(beta_arr)
f_arr = np.zeros_like(beta_arr)
eps_arr = np.zeros((len(shift_range),2))
tau_arr = np.zeros_like(eps_arr)
for i, p50 in enumerate(psi_T50_arr): 
    beta_arr[i], d_arr[i], f_arr[i], k = get_params(gX, p50, gS, psi_ssat, gC, psi_L90)
    eps_arr[i,0] = eps1(beta_arr[i], d_arr[i], f_arr[i], k)
    eps_arr[i,1] = eps2(beta_arr[i], d_arr[i], f_arr[i], k)
    tau_arr[i,0] = tau1(beta_arr[i], d_arr[i], f_arr[i], k)
    tau_arr[i,1] = tau2(beta_arr[i], d_arr[i], f_arr[i], k)

plt.figure()
plt.plot(beta_arr, eps_arr[:,0])
plt.plot(d_arr, eps_arr[:,0])
plt.plot(f_arr, eps_arr[:,0])
plt.figure()
plt.plot(beta_arr, tau_arr[:,0])
plt.plot(d_arr, tau_arr[:,0])
plt.plot(f_arr, tau_arr[:,0])
plt.show()

