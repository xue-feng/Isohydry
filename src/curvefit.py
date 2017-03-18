import xlrd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

filepath = '/Users/xuefeng/Dropbox/Projects/Isohydricity/hydraulic_traits.xls'
book = xlrd.open_workbook(filepath)
sheet = book.sheet_by_name('curvefit')
get_vals = lambda coln: np.asarray(filter(None, sheet.col_values(coln)[2:]), dtype='float')
label = ['vpd_pine', 'gmax_pine', 'vpd_juni', 'gmax_juni', 'psiL_juni', 'g_juni', 'psiL_pine', 'g_pine', 'psiX_juni', 'vc_juni', 'psiX_pine', 'vc_pine']
col_dict = {label[i]:get_vals(i) for i in range(12)}
# local variables dynamically assigned
locals().update(col_dict)

def exponential(x, a, b): return a*np.exp(-b*x)
def plc(P, P50, a): return 1.0 - 1.0/(1.0+np.exp(a*(P-P50)))
def expo_capped(x, b): return np.exp(-b*x)

''' for vpd dependencies ''' 
params_pine,_ = opt.curve_fit(exponential, vpd_pine, gmax_pine)
params_juni,_ = opt.curve_fit(exponential, vpd_juni, gmax_juni)
plt.figure()
plt.plot(vpd_pine, gmax_pine, 'o') 
plt.plot(vpd_juni, gmax_juni, 'o')
vpd_arr = np.arange(0,4,0.1)
plt.plot(vpd_arr, exponential(vpd_arr, params_pine[0], params_pine[1]))
plt.plot(vpd_arr, exponential(vpd_arr, params_juni[0], params_juni[1]))
print params_pine, params_juni

''' for vulnerability curves'''
params_pine,_ = opt.curve_fit(plc, psiX_pine, vc_pine)
params_juni,_ = opt.curve_fit(plc, psiX_juni, vc_juni)
plt.figure()
plt.plot(psiX_pine, vc_pine, 'o') 
plt.plot(psiX_juni, vc_juni, 'o')
psiX_arr = np.arange(0,12,0.1)
plt.plot(psiX_arr, plc(psiX_arr, params_pine[0], params_pine[1]))
plt.plot(psiX_arr, plc(psiX_arr, params_juni[0], params_juni[1]))
print params_pine, params_juni

''' for stomatal conductance ''' 
params_pine,_ = opt.curve_fit(expo_capped, psiL_pine, g_pine)
params_juni,_ = opt.curve_fit(expo_capped, psiL_juni, g_juni)
plt.figure()
plt.plot(psiL_pine, g_pine, 'o') 
plt.plot(psiL_juni, g_juni, 'o')
psiL_arr = np.arange(0,6,0.1)
plt.plot(psiL_arr, expo_capped(psiL_arr, params_pine))
plt.plot(psiL_arr, expo_capped(psiL_arr, params_juni))
print params_pine, params_juni
plt.show()