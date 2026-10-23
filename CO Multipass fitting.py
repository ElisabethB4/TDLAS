#%% 
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from pybaselines import Baseline, utils
import hapi as hp
from hapi import * 
from scipy.optimize import minimize, Bounds
from scipy.signal import find_peaks

def hapi_calculation(species,p, t, x, length, resltn, afwing, wnb, dnu):
            nu, coef = hp.absorptionCoefficient_Voigt(SourceTables=species,    # calculate voigt absorption coefficient
                                              WavenumberRange=(wnb),
                                              WavenumberStep=(dnu),
                                              Environment={'p':p, 'T':t}, 
                                              Diluent={'Air':(1-(x))}, 
                                              HITRAN_units=False)   
            coef_corrected = coef * x  # scale coefficent by mole fraction
            nu, absorp = hp.absorptionSpectrum(nu,coef_corrected,  # calculate absorption from abs coeff                 
                                              Environment={'l':length})  
            wl_convolve = 1e7/nu # calculate wavelength 
            wl_hapi, absorp_, i1, i2, slit = hp.convolveSpectrum(np.flip(wl_convolve),absorp, # convolve to insturment function   
                                              SlitFunction=hp.SLIT_GAUSSIAN,
                                              Resolution=resltn,AF_wing=afwing)
            wl = np.flip(wl_hapi) # wavelength
            
            return wl, absorp_ 

def objective(arguments, A_exp_wlc):
      
    # T, P, X, L = arguments[0], arguments[1], arguments[2], arguments[3] # pull t p x and l from the input array
    T = arguments[0] # pull t from the input array
    X = arguments[1] # pull t from the input array
    
    wl, A_sim = hapi_calculation('CO', P, T, X, length, resltn, afwing, wnb, dnu) # run Hcl calc for current condition

    A_sim_c = np.interp(wlc, np.flip(wl), np.flip(A_sim)) # interpolate to common wavelength scale`

    ssr = sum((A_exp_wlc - A_sim_c)** 2) # calculate r2

    return ssr
#%% 
cell_C = np.load('cell_C.npy')

abs1_C = np.load('abs1_C.npy')
abs1_5_C = np.load('abs1_5_C.npy')

abs2_C = np.load('abs2_C.npy')
abs2_5_C = np.load('abs2_5_C.npy')
x_narrow = np.load('x_narrow.npy')
#%%% 
hp.db_begin('HAPI_DATA') 

length = 2.54   # cm
resltn = 0.001 # cm-1
afwing = 0.0075  # instrument function width
scale = 0.7

wnb = [4290, 4293]
dnu = 0.0001 

P = 1
X = 1
T = 298

wl, cell_abs_sim = hapi_calculation('CO', P, T, X, length, resltn, afwing, wnb, dnu)

wl = np.flip(wl)
cell_abs_sim = np.flip(cell_abs_sim)
 
sort_idx = np.argsort(x_narrow)
xs = x_narrow[sort_idx]

wavelength_exp = 2329.8425 + (0.62) * (xs - xs.min()) / (xs.max() - xs.min()) # By hand guess
# wavelength_exp = 2329.8425 + (0.95) * (xs - xs.min()) / (xs.max() - xs.min())

wlc = np.arange(2329.95, 2330.5, 0.001)

A_exp_wlc = np.interp(wlc, wavelength_exp, cell_C)

A_sim_wlc = np.interp(wlc, wl, cell_abs_sim)

residuals = A_exp_wlc - A_sim_wlc
rsme_guess_1 = np.sqrt(np.sum((residuals ** 2) / len(residuals)))
print("Manual Cal RSME = {:.4f}".format(rsme_guess_1))

def wl_axis_cal(param, waveAxis, simAxis, abs_sim, abs_data):
    a1 = param[0]
 
    waveAxis_shift =  waveAxis + a1
    abs_sim_shift = np.interp(waveAxis_shift, simAxis, abs_sim)

    residuals = abs_data - abs_sim_shift
    residualrmse = np.sqrt(np.sum((residuals ** 2) / len(residuals)))  # root mean squared error
 
    return residualrmse

bnds = Bounds([-1e6], 
              [1e6])
init = [1e-3]

wl_fit_result = minimize(wl_axis_cal, init, bounds=bnds, tol= 1e-6, method='Nelder-Mead', args=(wavelength_exp, wl, cell_abs_sim, cell_C))
print("Fitted RMSE = ", wl_fit_result.fun)

wl_axis_fit = wavelength_exp + wl_fit_result.x[0]

A_exp_wlc = np.interp(wlc, wl_axis_fit, cell_C)
A_exp_wlc = A_exp_wlc 

plt.figure()
plt.plot(wlc, A_sim_wlc, label='Sim. T = {:.0f} K, X = {:.3f}'.format(T, X))
plt.plot(wlc, A_exp_wlc , label = 'Experimental')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Absorption')
# plt.xlim([2329.9,2330.5])
plt.title("Fitted Wl scale")
plt.legend()

#%% fit insturment function 
def fitfunct(arguments, A_exp_wlc):
      
    # T, P, X, L = arguments[0], arguments[1], arguments[2], arguments[3] # pull t p x and l from the input array
    resltn = arguments[0] # pull t from the input array
    afwing = arguments[1] *3 # pull t from the input array
    
    wl, A_sim = hapi_calculation('CO', P, T, X, length, resltn, afwing, wnb, dnu) # run Hcl calc for current condition

    A_sim_c = np.interp(wlc, np.flip(wl), np.flip(A_sim)) # interpolate to common wavelength scale`

    ssr = sum((A_exp_wlc - A_sim_c)** 2) # calculate r2

    return ssr

bnds = Bounds([-1e6, -1e6], 
              [1e6, 1e6])
init = [resltn, afwing]

fitfunct_result = minimize(fitfunct, init, bounds=bnds, tol= 1e-6, method='Nelder-Mead', args=(A_exp_wlc))

resltn = fitfunct_result.x[0] 
afwing = fitfunct_result.x[1] *3
#%% 

length = 2.54 * 7 # for multipass

abs_C = [abs1_C , abs1_5_C, abs2_C, abs2_5_C]
P = 1
eq_ratio =  [1,            1.5,          2,          2.5]
t        = [2364.88,      2186.52,      1898.70,    1643.98]
CO_x     = [0.01996,      0.12803,      0.20833,    0.27402] 

T_fit_array = []
X_fit_array = []
errorarray = []

abs_fitted_array = []
abs_CEA_array = []
abs_exp_array = []

for i in range(len(eq_ratio)):
    
    eq = eq_ratio[i] 
    T = t[i]
    X = CO_x[i]
    
    wl, absorp_sim = hapi_calculation('CO', P, T, X, length, resltn, afwing, wnb, dnu) # calculate absorp. for CEA values
    
    A_sim_wlc = np.interp(wlc, np.flip(wl), np.flip(absorp_sim))
    
    peaks_sim, _ = find_peaks(A_sim_wlc, distance=75, height= 0.005)
    
    spacing = np.diff(peaks_sim) - 2
    
    abs_CEA_array.append(A_sim_wlc)
    
    absorption_C = abs_C[i] 
    wl_axis = 2329.83 + (0.95 ) * (xs - xs.min()) / (xs.max() - xs.min())

    A_exp_wlc = np.interp(wlc, wl_axis, absorption_C) 
    
    bnds = Bounds([-1e6], 
                  [1e6])
    init = [1e-3]
    wl_fit_result = minimize(wl_axis_cal, init, bounds=bnds, tol= 1e-6, method='Nelder-Mead', args=(wl_axis, wlc, A_sim_wlc, absorption_C))
    wl_axis_fit =  wl_axis + wl_fit_result.x[0]

    A_exp_wlc = np.interp(wlc, wl_axis_fit, absorption_C) 
    abs_exp_array.append(A_exp_wlc)
    
    CEA_resid = sum((A_exp_wlc - A_sim_wlc)** 2) # calculate r2
    
    print(CEA_resid)

    bnds = Bounds([T-2000, 0], # lower bounds
                  [T+2000, 1]) # upper bounds
    
    initial_guess = [T, X ] # initial guess t=290K, p=1atm, x=0.3, l = 5 (given cell values)
    
    result = minimize(objective, initial_guess, bounds=bnds, tol=1e-4, method='Nelder-Mead', args=(A_exp_wlc)) # run minimizing routine
    
    print(result) # print minimizing routine results
    
    T_fit =  result.x[0] # read result for T
    X_fit = result.x[1] # read result for P
    
    T_fit_array.append(T_fit)
    X_fit_array.append(X_fit)
    error = result.fun
    errorarray.append(error)
    
    wl_min, A_minimized = hapi_calculation('CO', P, T_fit, X_fit, length, resltn, afwing, wnb, dnu)  # calculate absorb at minimized T,P,X,and L
    
    A_min_wlc = np.interp(wlc, np.flip(wl_min), np.flip(A_minimized)) # interpolate minimized to common scale
    
    abs_fitted_array.append(A_min_wlc)
    
    # plt.figure()
    # plt.plot(wlc, A_min_wlc, label='Fit T = {:.2f} K, X = {:.3f}'.format(T_fit, X_fit))
    # # plt.plot(wlc, A_sim_wlc, label='CEA T = {:.2f}, X = {:.2f}'.format(T, X))
    # plt.plot(wlc, A_exp_wlc , label='Exp.')
    # plt.xlabel('Wavelength (nm)')
    # plt.ylabel('Absorption Spectrum')
    # plt.title("Exp.and Fitted. Absorption Equivilance Ratio = {:.1f}".format(eq))
    # plt.legend()
    
abs_CEA_array = np.array(abs_CEA_array)
abs_exp_array = np.array(abs_exp_array)
abs_fitted_array = np.array(abs_fitted_array)
T_fit_array = np.array(T_fit_array)
X_fit_array = np.array(X_fit_array)
errorarray = np.array(errorarray)

np.savetxt('abs_exp_array.csv', (abs_exp_array), delimiter=',')
np.savetxt('abs_fitted_array.csv', (abs_fitted_array), delimiter=',')
np.savetxt('T_fit_array.csv', (T_fit_array), delimiter=',')
np.savetxt('X_fit_array.csv', (X_fit_array), delimiter=',')
np.savetxt('wavelength_array.csv', wlc, delimiter=',')


#%% 
for i in range(len(eq_ratio)):
    A_fitted_wlc = abs_fitted_array[i]
    A_exp_wlc = abs_exp_array[i]
    
    T_fit = T_fit_array[i]
    X_fit = X_fit_array[i]
    eq = eq_ratio[i]
    error = errorarray[i]
    
    plt.figure()
    plt.plot(wlc, A_fitted_wlc, label='Fit T = {:.2f} K, X = {:.2f}'.format(T_fit, X_fit))
    plt.plot(wlc, A_exp_wlc , label='Exp.')
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Absorption Spectrum')
    plt.title("Exp., CEA guess, and Fitted. Absorption Equivilance Ratio = {:.1f}".format(eq))
    plt.legend()
    plt.rcParams["text.usetex"] = False

    nu = 1e7 / wlc  
    fit_transmission = 1 - A_fitted_wlc
    data_transmission = 1 - A_exp_wlc

    fig, axs = plt.subplots(2, 1, figsize=[10, 6])
    axs[0].plot(
        nu,
        fit_transmission,
        color="tab:blue",
        label=rf"Spectra at $X = {X_fit:.4f},\ T = {T_fit:.1f}\ \mathrm{{K}},\ L = {length:.1f}\ \mathrm{{cm}}$",
    )
    axs[0].plot(
        nu[::3],
        data_transmission[::3],
        marker="o",
        linestyle="none",
        color="red",
        label ="Data",
    )

    axs[0].set_ylabel(r"Transmission")
    axs[0].legend()
    axs[0].grid(True)
    axs[0].set_ylim(0.9, 1.01)
    axs[0].set_xlabel(r"Wavenumber, $(\mathrm{cm}^{-1})$")
    axs[0].ticklabel_format(axis="x", style="plain", useOffset=False)

    axs[1].plot(nu, data_transmission - fit_transmission, color="tab:blue", label=f"$RSME = {error**2:.2E}$")
    axs[1].set_ylabel(r"Residuals")
    axs[1].set_xlabel(r"Wavenumber, $(\mathrm{cm}^{-1})$")
    axs[1].grid(True)
    axs[1].legend(loc="lower left")
    axs[1].ticklabel_format(axis="x", style="plain", useOffset=False)
    fig.tight_layout(h_pad=2.0)

    

#%%

plt.rcParams["text.usetex"] = False

fig_summary, axs = plt.subplots(2, 1)
axs[0].plot(eq_ratio, T_fit_array, marker="o", linestyle="none", color="red")
axs[0].set_ylabel("Fitted temperature (K)")

axs[1].plot(eq_ratio, X_fit_array, marker="o", linestyle="none", color="blue")
axs[1].set_ylabel("Fitted CO mole fraction")
axs[1].set_xlabel("Equivalence ratio")

fig_summary.tight_layout(h_pad=2.0)
plt.show()

# %%
