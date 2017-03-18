import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
import scipy.optimize as opt
import scipy.interpolate as interpolate
import scipy.integrate as integ
import time

def params(lam, gam, eta, k, sw, sst, s1):
    '''need to make sure R!=0'''
    R = gam*eta-lam
    pi = lam/eta*(s1-sst)
    theta = lam/eta*(sst-sw)
    zeta = gam*eta/k*(1-s1)
    chi = lam/k*(1-s1)
    G1 = sp.gammaincc(chi, zeta)*sp.gamma(chi)
    G2 = sp.gammaincc(chi, zeta+gam*(1-s1))*sp.gamma(chi)
    G3 = sp.gammainc(theta, gam*(sst-sw))*sp.gamma(theta)
    G4 = sp.gammaincc(chi+1, zeta)*sp.gamma(chi+1)
    G5 = sp.gammaincc(chi+1, zeta+gam*(1-s1))*sp.gamma(chi+1)
    G6 = sp.gammainc(theta+1, gam*(sst-sw))*sp.gamma(theta+1)
    return R, pi, theta, zeta, chi, G1, G2, G3, G4, G5, G6

def normC(lam, gam, eta, k, sw, sst, s1):
    R, pi, theta, zeta, chi, G1, G2, G3, _, _, _ = params(lam, gam, eta, k, sw, sst, s1)
    C = 1.0/( (np.exp(-gam*sst)-np.exp(-gam*s1+pi))/R 
          + 1.0/eta* np.exp(-gam*sw)* (sst-sw)**(-theta+1)* gam**(-theta)* G3
          - 1.0/k* np.exp(-gam*s1+zeta+pi)* (1-s1)* zeta**(-chi)* (G2-G1)
         )
    return C

def smu(lam, gam, eta, k, sw, sst, s1):
    R, pi, theta, zeta, chi, G1, G2, G3, G4, G5, G6 = params(lam, gam, eta, k, sw, sst, s1)
    C = normC(lam, gam, eta, k, sw, sst, s1)
    s_mean = C*( (np.exp(-gam*sst)*(eta+sst*R)-np.exp(-gam*s1+pi)*(eta+s1*R))/R**2
                 + 1.0/k* np.exp(-gam*s1+zeta+pi)*(1-s1)*zeta**(-chi)* ( (s1-eta/k*(1-s1))*(G1-G2)+ 1.0/gam*(G4-G5))
                 + 1/eta* np.exp(-gam*sw)* (sst-sw)**(-theta+1)* gam**(-theta)* (sw*G3 + G6/gam)
                ) 
    return s_mean

def ps_arr(s, lam, gam, eta, k, sw, sst, s1):
    ''' outputs values for p(s) at steady state given input parameters  '''
    p_s = np.zeros_like(s)
    C = normC(lam, gam, eta, k, sw, sst, s1)
    ilow = (s>sw)&(s<=sst)
    imid = (s>sst)&(s<=s1)
    ihigh = (s>s1)&(s<=1.0)
    p_s[ilow] = C/eta*np.exp(-gam*s[ilow])* ((s[ilow]-sw)/(sst-sw))**(lam/eta*(sst-sw)-1)
    p_s[imid] = C/eta*np.exp(-gam*s[imid])* np.exp(lam/eta*(s[imid]-sst))
    p_s[ihigh] = C/eta*np.exp(-gam*s[ihigh])* np.exp(lam/eta*(s1-sst))* (k/eta*(s[ihigh]-s1)/(1-s1)+1)**(lam/k*(1-s1)-1)
    return p_s

def ps(s, lam, gam, eta, k, sw, sst, s1):
    ''' p(s) based on single value of s '''
    C = normC(lam, gam, eta, k, sw, sst, s1)
    if (s>sw)&(s<=sst): 
        p_s = C/eta*np.exp(-gam*s)* ((s -sw)/(sst-sw))**(lam/eta*(sst-sw)-1)
    ilow = (s>sw)&(s<=sst)
    imid = (s>sst)&(s<=s1)
    ihigh = (s>s1)&(s<=1.0)
    p_s[ilow] = C/eta*np.exp(-gam*s[ilow])* ((s[ilow]-sw)/(sst-sw))**(lam/eta*(sst-sw)-1)
    p_s[imid] = C/eta*np.exp(-gam*s[imid])* np.exp(lam/eta*(s[imid]-sst))
    p_s[ihigh] = C/eta*np.exp(-gam*s[ihigh])* np.exp(lam/eta*(s1-sst))* (k/eta*(s[ihigh]-s1)/(1-s1)+1)**(lam/k*(1-s1)-1)
    return p_s

def rho(s, lam, gam, eta, k, sw, sst, s1):
    if s<=sw: rho_s = 0.0 
    elif (s>sw)&(s<=sst): rho_s = eta*(s-sw)/(sst-sw)
    elif (s>sst)&(s<=s1): rho_s = eta
    elif (s>s1)&(s<=1.0): rho_s = eta + k*(s-s1)/(1.0-s1)
    return rho_s

def gcorr(s, lam, gam, eta, k, sw, sst, s1, Zr,corr_params):
    c0,c1,c2,Zr0 = corr_params
    if s<=sw: eps = 0 
    elif (s>sw)&(s<=sst): eps = c0 + c1*(s-sw)/(sst-sw)
    elif (s>sst)&(s<=s1): eps = c0 + c1 + c2*np.sqrt(Zr0/Zr)*(s-sst)/(s1-sst)
    elif (s>s1)&(s<=1.0): eps = c0 + c1 + c2*np.sqrt(Zr0/Zr)
    return lam/gam * eps

def rho_gcorrected(s, lam, gam, eta, k, sw, sst, s1, Zr, corr_params):
    c0,c1,c2,Zr0 = corr_params
    if s<=sw: eps = 0.0; rho_s = 0.0
    elif (s>sw)&(s<=sst): eps = c0 + c1*(s-sw)/(sst-sw);                       rho_s = eta*(s-sw)/(sst-sw) 
    elif (s>sst)&(s<=s1): eps = c0 + c1 + c2*np.sqrt(Zr0/Zr)*(s-sst)/(s1-sst); rho_s = eta 
    elif (s>s1)&(s<=1.0): eps = c0 + c1 + c2*np.sqrt(Zr0/Zr);                  rho_s = eta + k*(s-s1)/(1.0-s1)
    return rho_s + lam/gam * eps

# need inverse transform sampling to generate samples according to p(s)
def inverse_transform_sampling(func, params, n_samples=1000):
    sunif = np.linspace(1e-10, 1.0+1e-10, n_samples)
    pdf_vals = func(sunif, *params)
    dx = 1.0/len(sunif)
    cdf_vals = dx*np.cumsum(pdf_vals)
    inv_cdf = interpolate.interp1d(cdf_vals, sunif, bounds_error=False, fill_value=(sw,1.0))
    r = np.random.rand(n_samples)
    return inv_cdf(r)

def smu_inverse_func(lam, gam, eta, k, sw, sst, s1, s_mean):
    # helper function for finding inverse of smu as function of lam
    # first parameter used to invert
    return np.abs(s_mean - smu(lam, gam, eta, k, sw, sst, s1))

def calculate_matching_lam(s_mean,gam,eta,k,sw,sst,s1):
    lam_opt = opt.fminbound(smu_inverse_func, 0.0, 1.0, xtol = 1e-08, args=(gam, eta, k, sw, sst, s1, s_mean))
    return lam_opt

def dsmu_dt_Laio2002(s_mean, t, lam, gam, eta, k, sw, sst, s1, Zr, corr_params):
    dsdt = lam/gam - rho_gcorrected(s_mean, lam, gam, eta, k, sw, sst, s1, Zr, corr_params)
    return dsdt

def dsmu_dt_uncorrected(s_mean, t, lam, gam, eta, k, sw, sst, s1, Zr, corr_params):
    dsdt = lam/gam - rho(s_mean, lam, gam, eta, k, sw, sst, s1)
    return dsdt
 
def dsmu_dt_selfconsistent(s_mean, t, lam, gam, eta, k, sw, sst, s1, Zr, corr_params):
    ''' input s_mean need to be a scalar because inverse can be performed only on a scalar '''
    # find the inverse of smu (e.g., the lambda parameter giving the desired s)
    lam_opt = calculate_matching_lam(s_mean,gam,eta,k,sw,sst,s1)
    # calculate s based on lam 
    s_sampled = inverse_transform_sampling(ps_arr,(lam_opt, gam, eta, k, sw, sst, s1), 1000)
    Qmu = lam/gam * np.mean( np.exp(-gam*(1.0-s_sampled)) )
    rho_corrected = rho_gcorrected(s_mean, lam, gam, eta, k, sw, sst, s1, Zr, corr_params)
    dsdt = lam/gam - rho_corrected - Qmu
    return dsdt

def runge_kutta4( f, tmax, nsteps, IV, params ):
    h = (tmax)/float(nsteps)                 # determine step-size
    t = np.arange( 0.0, tmax+h, h )          # create mesh
    w = np.zeros((nsteps+1,))                # initialize w
    t[0], w[0] = IV                     # set initial values
    for i in range(1,nsteps+1):              # apply Fourth Order Runge-Kutta Method
        k1 = h * f( w[i-1], t[i-1], *params )
        k2 = h * f( w[i-1] + k1 / 2.0, t[i-1] + h / 2.0, *params)
        k3 = h * f( w[i-1] + k2 / 2.0, t[i-1] + h / 2.0, *params)
        k4 = h * f( w[i-1] + k3, t[i], *params )
        w[i] = w[i-1] + ( k1 + 2.0 * k2 + 2.0 * k3 + k4 ) / 6.0
    return w

# s_arr = np.arange(0.6, 1.0, 0.01)
# lam_params = np.zeros_like(s_arr)
# smu_arr = np.zeros_like(s_arr)
# for i, s in enumerate(s_arr):
#     lam_params[i] = calculate_matching_lam(s, gam, eta, k, sw, sst, s1)
#     smu_arr[i] = np.mean( inverse_transform_sampling(ps_arr,(lam_params[i], gam, eta, k, sw, sst, s1), 10000) )
# plt.plot(s_arr, lam_params)
# plt.plot(s_arr, smu_arr)
# plt.show()
# 
# lam_params = np.arange(0.1,0.5,0.05)
# for j,l in enumerate(lam_params): 
#     plt.figure()
#     plt.hist(inverse_transform_sampling(ps_arr,(l, gam, eta, k, sw, sst, s1), 1000))
# plt.show()

# st_Laio2002 = integ.odeint(dsmu_dt_Laio2002,xInit,tRun,args=(lam,gam,eta,k,sw,sst,s1,Zr)).T[0]
# st = integ.odeint(dsmu_dt,xInit,tRun,args=(lam,gam,eta,k,sw,sst,s1,Zr)).T[0]
# plt.plot(tRun, st)
# plt.plot(tRun, st_Laio2002)
# plt.show()

''' what's left: 
1. simulate time dependent soil moisture pdfs to compare with self-consistent solution.
2. look into why the solution isn't well defined for Zr<25.0cm ??  
'''

def simulate_s_t(depths, tRun, sInit, lam, gam, eta, k, sw, sst, s1):
    ''' simulate for a single trajectory'''
    s_t = np.zeros(len(tRun))
    s0 = sInit
    for i in range(len(tRun)): 
        R_normed = depths[i]
        Infil_normed = min(R_normed, 1.0-s0)
        ET_L_normed = rho(s0, lam, gam, eta, k, sw, sst, s1) # deterministic
        s_out = s0 + Infil_normed - dt*(ET_L_normed)
        s_t[i] = s_out
        s0 = s_out
    return s_t

def simulate_ps_t(n_trajectories, tRun, s0, lam, gam, eta, k, sw, sst, s1):
    ''' simulate for multiple trajectories over time '''
    # generating rainfall depth all at once
    size = len(tRun)*n_trajectories
    depthExp = -np.log(1.0-np.random.random(size=size))/gam
    freqUnif = np.random.random(size=size)
    
    depth = np.zeros(size)
    # the occurence of rainfall in any independent interval is lam*dt
    depth[freqUnif<np.tile(lam,size)*dt] = depthExp[freqUnif<np.tile(lam,size)*dt] # rain falls according to prob within an increment
    depth_re = np.reshape(depth, (n_trajectories, len(tRun)))
    ps_samples = np.zeros((n_trajectories, len(tRun)))
    for nsim in range(n_trajectories):
        s_t = simulate_s_t(depth_re[nsim],tRun,s0,lam, gam, eta, k, sw, sst, s1)
        ps_samples[nsim] = s_t
    return ps_samples

## plot trajectory of mean s ## 
sInit=1.0; dt = 0.1; tmax=50
tRun = np.arange(0,tmax+dt,dt); nsteps=int(tmax/dt)
s0=sInit; n_trajectories = 1000

sunif = np.arange(0,1.0,0.001)
lam_arr = np.linspace(0.01,1.0,1000); lam=0.4
eta = 0.05; gam = 25; k = 0.1; sw=0.11; sst = 0.31; s1=0.6; Zr=50.0 # Zr in cm

c0=0; c1= -0.2; c2=0.75; Zr0=50 # As in Laio2002
c0=0; c1= -0.3; c2=0.5; Zr0=35 # As in Laio2002, works better when s<sst 
corr_params=(c0,c1,c2,Zr0) 
int_params = (lam, gam, eta, k, sw, sst, s1, Zr, corr_params)

print 'start simulating'
ps_samples = simulate_ps_t(n_trajectories,tRun,s0,lam, gam, eta, k, sw, sst, s1)

# lam = 0.2; Zr=50.0; gam=5.5
# ps_samples = simulate_ps_t(n_trajectories,tRun,s0,lam, gam, eta, k, sw, sst, s1)
# test_params = lambda corr_params: integ.odeint(dsmu_dt_Laio2002,sInit,tRun,args=(lam,gam,eta,k,sw,sst,s1,Zr,corr_params)).T[0]
# param_list = [(0,-0.2,0.75,50), (0,-0.2,0.5,35)]
# for cp in param_list:
#     plt.plot(tRun, test_params(cp))
# plt.plot(tRun, np.mean(ps_samples,axis=0), lw=2, color='black')
# plt.show()    

print 'start modeling'
start_time = time.time()
st_selfconsist = runge_kutta4(dsmu_dt_selfconsistent, tmax, nsteps, (0.0, sInit), int_params)
st_Laio2002 = integ.odeint(dsmu_dt_Laio2002,sInit,tRun,args=(lam,gam,eta,k,sw,sst,s1,Zr,corr_params)).T[0]
print (time.time() - start_time)/60.0 

print 'start plotting'
for nsim in range(20):
    plt.plot(tRun, ps_samples[nsim], lw=0.5, color='lightgray')
plt.plot(tRun, ps_samples[0], lw=0.5, color='black')
plt.plot(tRun, np.mean(ps_samples,axis=0), lw=2, color='black')
plt.plot(tRun, st_selfconsist, lw=2, color='blue')
plt.plot(tRun, st_Laio2002, lw=2, color='green')
plt.plot(tRun, np.zeros(len(tRun))*(sw + lam/(gam*eta)*(sst-sw)) )
plt.show()


''' why does the decay not go down to steady state value???
This is because lam/gam is "buffering" the decay. 
Need to approximate what s_steady_state is, then adjust the correction factor to match that, 
assuming once the trajectory has decayed to that point, leakage plays very little role 

sterm = sw + lam/(gam*eta)*(sst-sw) --> where does this come from??

sw*(1.0 - lam/(gam*eta)) + lam/(gam*eta)*sst

dsdt = lam/gam - eta
'''
# f_corr = np.zeros(len(sunif))
# rho_arr = np.zeros_like(f_corr)
# for i,s in enumerate(sunif):
#     f_corr[i] = rho_gcorrected(s, lam, gam, eta, k, sw, sst, s1, Zr, corr_params)
#     rho_arr[i] = rho(s, lam, gam, eta, k, sw, sst, s1)
# plt.figure()
# plt.plot(sunif, f_corr,'--')
# plt.plot(sunif, rho_arr)
# plt.show()