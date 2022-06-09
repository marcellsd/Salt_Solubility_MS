# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 14:23:12 2022

@author: marcell
"""
import sys

sys.path.append('Modules')

import pandas as pd
from numpy import array, linspace, log, empty, exp, mean, concatenate, atleast_2d
import matplotlib.pyplot as plt
from tkinter import filedialog, Tk

from Modules import transformations
from Modules import pitzer
from Modules import excess
from Modules import statistics
from matplotlib.backends.backend_pdf import PdfPages



root = Tk()
root.withdraw()
data_file_name = filedialog.askopenfilename()

df = pd.read_excel(data_file_name, index_col=('Parameter'))
Rm = df.loc['Cation radius', 'Value']
Rx = df.loc['Anion radius', 'Value']
Zm = df.loc['Cation charge', 'Value']
Zx = df.loc['Anion Charge', 'Value']
nu_M = df.loc['Cation stoichiometric coefficient', 'Value']
nu_X = df.loc['Anion stoichiometric coefficient', 'Value']
mm_MX = df.loc['Salt Molar Mass (g/mol)', 'Value']
MM_MX = mm_MX/1000 #Kg/mol
DeltaG_transf = df.loc['Transf DeltaG (J/mol)', 'Value']
cation_name = df.loc['Cation', 'Value']
anion_name = df.loc['Anion', 'Value']

if nu_M==1 and nu_X==1:
    salt_name = cation_name + anion_name
elif nu_M==1 and nu_X!=1:
    salt_name = cation_name + anion_name + str(nu_X)
elif nu_M!=1 and nu_X==1:
    salt_name = cation_name + str(nu_M) + anion_name
else:
    salt_name = cation_name + str(nu_M) + anion_name + str(nu_X)

#Pitzer Coefficients
Beta_0 = df.loc['Beta0_ref_25C', 'Value']
Beta_1 = df.loc['Beta1_ref_25C', 'Value']
Cphi = df.loc['Cphi_ref_25C', 'Value']

data_file_sol = pd.read_excel(data_file_name, 'Data')
T = pd.DataFrame(data_file_sol, columns=['T (K)']).to_numpy()
w_H2O_SF = pd.DataFrame(data_file_sol, columns=['w_H20_SF']).to_numpy()
w_MX_all = pd.DataFrame(data_file_sol, columns=['w_salt']).to_numpy()
T_filter = pd.DataFrame(data_file_sol, columns=['T (K)']).drop_duplicates().to_numpy()

#if w_H2O_SF ≃ 1 -> Water-Salt binary
#if w_H2O_SF ≃ 0 -> Salt-Meg binary

#Separation of data into ternary and binary
w_MX_H2O_bin_L, w_H2O_bin_L = [], []
w_MEG_bin_L, w_MX_MEG_bin_L = [], []
w_H2O_tern_L, w_MEG_tern_L, w_MX_tern_L, w_H2O_SF_tern_L, w_MEG_SF_tern_L = [], [], [], [], []
w_H2O_All_L, w_MEG_All_L, w_MX_All_L, w_H2O_SF_All_L, w_MEG_SF_All_L, T_All_L = [], [], [], [], [], []   
w_H2O_SF_in_bin_H2O_L, w_H2O_SF_in_bin_MEG_L = [], []
w_MEG_SF_in_bin_MEG_L, w_MEG_SF_in_bin_H2O_L = [], []


vector_length = len(w_H2O_SF)
for i in range (0, vector_length):
    w_H2O_All_L.append(w_H2O_SF[i]*(1 - w_MX_all[i]))
    w_MX_All_L.append(w_MX_all[i])
    w_MEG_All_L.append(1 - w_H2O_SF[i]*(1-w_MX_all[i]) - w_MX_all[i])
    w_H2O_SF_All_L.append(w_H2O_SF[i])
    w_MEG_SF_All_L.append(1-w_H2O_SF[i])
    T_All_L.append(T[i])
    
    if w_H2O_SF[i] >= 0.99:
        w_MX_H2O_bin_L.append(w_MX_all[i])
        w_H2O_bin_L.append(w_H2O_SF[i]*(1 - w_MX_all[i]))
        w_H2O_SF_in_bin_H2O_L.append(w_H2O_SF[i])
        w_MEG_SF_in_bin_H2O_L.append(1-w_H2O_SF[i])
    
    elif w_H2O_SF[i] <= 0.01:
        w_MEG_bin_L.append((1-w_H2O_SF[i])*(1 - w_MX_all[i]))
        w_MX_MEG_bin_L.append(w_MX_all[i])
        w_MEG_SF_in_bin_MEG_L.append(1-w_H2O_SF[i])
        w_H2O_SF_in_bin_MEG_L.append(w_H2O_SF[i])
    else:
        w_H2O_tern_L.append(w_H2O_SF[i]*(1 - w_MX_all[i]))
        w_MX_tern_L.append(w_MX_all[i])
        w_MEG_tern_L.append(1 - w_H2O_SF[i]*(1-w_MX_all[i]) - w_MX_all[i])
        w_H2O_SF_tern_L.append(w_H2O_SF[i])
        w_MEG_SF_tern_L.append(1-w_H2O_SF[i])  

#Substitute values of zero

for i in range (0, len(w_MX_MEG_bin_L)):
    if w_MX_MEG_bin_L[i] == 0:
        w_MX_MEG_bin_L[i] = 1e-6

#List to array
w_MX_H2O_bin, w_H2O_bin  = array(w_MX_H2O_bin_L), array(w_H2O_bin_L)
w_MEG_bin, w_MX_MEG_bin = array(w_MEG_bin_L, dtype="object"), array(w_MX_MEG_bin_L, dtype="object")
w_H2O_tern, w_MEG_tern, w_MX_tern, w_H2O_SF_tern = array(w_H2O_tern_L), array(w_MEG_tern_L), array(w_MX_tern_L), array(w_H2O_SF_tern_L)
w_H2O_All, w_MEG_All, w_MX_All, w_H2O_SF_All, w_MEG_SF_All, T_All = array(w_H2O_All_L), array(w_MEG_All_L), array(w_MX_All_L), array(w_H2O_SF_All_L), array(w_MEG_SF_All_L), array(T_All_L)
w_MEG_SF_tern = array(w_MEG_SF_tern_L)
w_H2O_SF_in_bin_H2O, w_MEG_SF_in_bin_H2O, w_MEG_SF_in_bin_MEG, w_H2O_SF_in_bin_MEG = array(w_H2O_SF_in_bin_H2O_L), array(w_MEG_SF_in_bin_H2O_L), array(w_MEG_SF_in_bin_MEG_L), array(w_H2O_SF_in_bin_MEG_L)


#Mean Calculations
w_MX_H2O_bin_mean = array(mean(w_MX_H2O_bin))
w_MX_MEG_bin_mean = array(mean(w_MX_MEG_bin))
w_H2O_SF_in_bin_H2O_mean = array(mean(w_H2O_SF_in_bin_H2O))
w_MEG_SF_in_bin_H2O_mean = array(mean(w_MEG_SF_in_bin_H2O))
w_MEG_SF_in_bin_MEG_mean = array(mean(w_MEG_SF_in_bin_MEG))
w_H2O_SF_in_bin_MEG_mean = array(mean(w_H2O_SF_in_bin_MEG))


#Concatenating Arrays
w_H2O_SF_All_True = concatenate([atleast_2d(a) for a in [w_H2O_SF_in_bin_H2O_mean,w_H2O_SF_tern,w_H2O_SF_in_bin_MEG_mean]])
w_MEG_SF_All_True = concatenate([atleast_2d(a) for a in [w_MEG_SF_in_bin_H2O_mean,w_MEG_SF_tern,w_MEG_SF_in_bin_MEG_mean]])
w_MX_All_True = concatenate([atleast_2d(a) for a in [w_MX_H2O_bin_mean,w_MX_tern,w_MX_MEG_bin_mean]])

for i in range (0, len(w_H2O_SF_All_True)):
    if w_H2O_SF_All_True[i] >= 1e-3:
        w_H2O_SF_All_True[i] == 0
for i in range (0, len(w_MEG_SF_All_True)):
    if w_MEG_SF_All_True[i] >= 1e-3:
        w_MEG_SF_All_True[i] == 0


#Transforming Mass fraction into Molar Salt-Free

x_MEG_SF_All = transformations.x_MEG_SF(w_H2O_SF_All_True)
x_H2O_SF_All = (1-x_MEG_SF_All)

x_MEG_SF_tern = transformations.x_MEG_SF(w_H2O_SF_tern)
x_H2O_SF_tern = (1-x_MEG_SF_tern)


#Molalities Calculations

b_MX_H2O_bin = transformations.mass_fraction_to_molality(MM_MX, w_MX_H2O_bin)
b_MX_MEG_bin = transformations.mass_fraction_to_molality(MM_MX, w_MX_MEG_bin)
b_MX_H2O_mean = array(mean(b_MX_H2O_bin))
b_MX_MEG_mean = array(mean(b_MX_MEG_bin))


#Pitzer

A_phi_H2O = pitzer.A_phi_H2O(T_filter)

rho_H2O, rho_MEG, epsilon_r_H2O, epsilon_r_MEG = pitzer.rho_and_epsilon_r_of_pures(T_filter)

A_phi_MEG = pitzer.A_phi_MEG(A_phi_H2O, rho_H2O, epsilon_r_H2O, rho_MEG, epsilon_r_MEG)

rho_H2O_MEG, epsilon_r_H2O_MEG = pitzer.rho_and_epsilon_r_of_mixing(rho_H2O, rho_MEG, epsilon_r_H2O, epsilon_r_MEG, 
                                                                    w_H2O_SF_All_True, x_MEG_SF_All, x_H2O_SF_All, w_MEG_SF_All_True, T_filter)
A_phi_H2O_MEG = pitzer.A_phi_of_Mixing(rho_H2O_MEG, rho_H2O, epsilon_r_H2O_MEG, epsilon_r_H2O, T_filter)

Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O = pitzer.pitzerParameters_MX_H2O(Beta_0, Beta_1, Cphi, T_filter, Zm, Rm, Rx)

ln_gamma_MX_H2O = pitzer.ln_gamma_MX_SI(nu_M, nu_X, Zm, Zx, b_MX_H2O_mean, A_phi_H2O, Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O)

Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG = pitzer.pitzerParameter_MX_MEG_SI(nu_M, nu_X, A_phi_MEG, Beta_0_MX_H2O, Beta_1_MX_H2O,
                                                                          Cphi_MX_H2O, ln_gamma_MX_H2O, b_MX_H2O_mean, b_MX_MEG_mean, 
                                                                          T_filter, Zm, Zx, DeltaG_transf)
Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, Cphi_MX_H2O_MEG = pitzer.Pitzer_parameters_in_mixtures_SI(Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O, 
                                                                                             Beta_1_MX_MEG, epsilon_r_H2O, epsilon_r_MEG, epsilon_r_H2O_MEG)


# Excess chemical potential calculations

b_MX_H2O_MEG = transformations.mass_fraction_to_molality(MM_MX, w_MX_All_True)
ln_gamma_MX_MEG = pitzer.ln_gamma_MX_SI(nu_M, nu_X, Zm, Zx, b_MX_MEG_mean, A_phi_MEG, Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG)
ln_gamma_MX_H2O_MEG = pitzer.ln_gamma_MX_SI(nu_M, nu_X, Zm, Zx, b_MX_H2O_MEG, A_phi_H2O_MEG, Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, Cphi_MX_H2O_MEG, False)  
ln_gamma_MX_H2O_MEG_ideal = excess.ln_gamma_MX_H2O_MEG_ideal(x_MEG_SF_All, ln_gamma_MX_H2O, ln_gamma_MX_MEG)
ln_gamma_MX_H2O_MEG_exc = excess.ln_gamma_MX_H2O_MEG_exc(ln_gamma_MX_H2O_MEG_ideal, ln_gamma_MX_H2O_MEG)
ln_b_MX_H2O_MEG = log(b_MX_H2O_MEG)
ln_b_MX_H2O_MEG_ideal = excess.ln_b_MX_H2O_MEG_ideal(x_MEG_SF_All ,b_MX_H2O_mean, b_MX_MEG_mean)
ln_b_MX_H2O_MEG_exc = excess.ln_b_MX_H2O_MEG_exc(ln_b_MX_H2O_MEG_ideal, ln_b_MX_H2O_MEG)
exc_chemical_potential = excess.excess_chemical_potential(ln_b_MX_H2O_MEG_exc, ln_gamma_MX_H2O_MEG_exc)

#Number of parameters in chemical Potential model (choose 3 ou 2)
n = 3
     
Thetas_pot = excess.n_linear_regression_chemical_potential_SI(exc_chemical_potential, x_MEG_SF_All, n)
pot_calc = excess.model_excess_chemical_potential_SI(x_MEG_SF_All, Thetas_pot, n)
x_MEG = linspace( 0.0,1.0, 100, endpoint=True)
pot_calc_plot = excess.model_excess_chemical_potential_SI(x_MEG, Thetas_pot, n)


#Chemical Potential statistics
R_2_pot = statistics.coefficient_of_determination2(exc_chemical_potential, pot_calc)
AAD_pot, maxAD_pot, AARD_pot, maxARD_pot = statistics.DAM_DMR(exc_chemical_potential, pot_calc)


#Mixture solubility calculations
number_of_points = 31

x_MEG_SF_opt = linspace(0.0,1.0,number_of_points,endpoint=True)
x_H2O_SF_opt = 1-x_MEG_SF_opt
w_MEG_SF_opt = transformations.w_MEG_SF(x_MEG_SF_opt)
w_H2O_SF_opt = 1-w_MEG_SF_opt
b_MX_H2O_MEG_optimizated = empty(number_of_points)
gamma_MX_H2O_MEG_optimizated = empty(number_of_points)

b_MX_H2O_MEG_optimizated, gamma_MX_H2O_MEG_optimizated = excess.optimization_b_MX_H2O_MEG_SI(x_MEG_SF_opt, w_H2O_SF_opt, x_H2O_SF_opt, w_MEG_SF_opt, T_filter[0], 
                                                                                                     Thetas_pot, b_MX_H2O_mean, b_MX_MEG_mean, 
                                                                                                     Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O,
                                                                                                     Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG,
                                                                                                     Zm, nu_M, nu_X, Zx, n)

b_MX_H2O_MEG_est = empty(number_of_points)
gamma_MX_H2O_MEG_est = empty(number_of_points)
b_MX_H2O_MEG_est, gamma_MX_H2O_MEG_est = excess.optimization_b_MX_H2O_MEG_SI(x_MEG_SF_All[:,0], w_H2O_SF_All_True[:,0], x_H2O_SF_All[:,0], w_MEG_SF_All_True[:,0], T_filter[0], 
                                                                                                     Thetas_pot, b_MX_H2O_mean, b_MX_MEG_mean, 
                                                                                                     Beta_0_MX_H2O[0], Beta_1_MX_H2O[0], Cphi_MX_H2O[0],
                                                                                                     Beta_0_MX_MEG[0], Beta_1_MX_MEG, Cphi_MX_MEG[0],
                                                                                                     Zm, nu_M, nu_X, Zx, n)
gamma_MX_H2O_MEG = exp(ln_gamma_MX_H2O_MEG)


#Solubility and Gamma Statistics
R_2_b=statistics.coefficient_of_determination2(b_MX_H2O_MEG, b_MX_H2O_MEG_est)
R_2_gamma=statistics.coefficient_of_determination2(gamma_MX_H2O_MEG, gamma_MX_H2O_MEG_est)
AAD_b, maxAD_b, AARD_b, maxARD_b = statistics.DAM_DMR(b_MX_H2O_MEG, b_MX_H2O_MEG_est)
AAD_g, maxAD_g, AARD_g, maxARD_g = statistics.DAM_DMR(gamma_MX_H2O_MEG, gamma_MX_H2O_MEG_est)

print (f"R_2_b: {R_2_b:.4f}, R_2_gamma: {R_2_gamma:.4f}, R_Pot: {R_2_pot:.4f}")

fig1 = plt.figure(1)      
plt.figure(1)
plt.xlabel("x_MEG_SF")
plt.ylabel("Potential")
plt.plot(x_MEG_SF_All, exc_chemical_potential, '.')
plt.plot(x_MEG, pot_calc_plot, label=T_filter)
plt.title("Calculated Chemical Potentials")
plt.legend()


#Graphic b
fig2 = plt.figure(2)
plt.xlabel("x_MEG_SF")
plt.ylabel("b")
plt.plot(x_MEG_SF_All, b_MX_H2O_MEG, '.')
plt.plot(x_MEG_SF_opt, b_MX_H2O_MEG_optimizated, label=T_filter)
plt.title("Calculated b_MX_H2O_MEG")
plt.legend()

#Graphic gamma
fig3 = plt.figure(3)
plt.xlabel("x_MEG_SF")
plt.ylabel("gamma")
plt.plot(x_MEG_SF_All,gamma_MX_H2O_MEG, '.')
plt.plot(x_MEG_SF_opt, gamma_MX_H2O_MEG_optimizated, label=T_filter)
plt.title("Calculated gamma_MX_H2O_MEG")
plt.legend()

fig4 = plt.figure(4)
plt.xlabel("x_MEG_SF")
plt.ylabel("ln_gamma")
plt.plot(x_MEG_SF_All,ln_gamma_MX_H2O_MEG, '.')
plt.plot(x_MEG_SF_opt, log(gamma_MX_H2O_MEG_optimizated), label=T_filter)
plt.title("Calculated gamma_MX_H2O_MEG")
plt.legend()

fig5 = plt.figure(5)
plt.xlabel("x_MEG_SF")
plt.ylabel("ln_b")
plt.plot(x_MEG_SF_All, ln_b_MX_H2O_MEG, '.')
plt.plot(x_MEG_SF_opt, log(b_MX_H2O_MEG_optimizated), label=T_filter)
plt.title("Calculated b_MX_H2O_MEG")
plt.legend()

pp = PdfPages(f'output_{salt_name}_{int(T_filter)}_SI.pdf')
plt.savefig(pp, format='pdf')
pp.savefig(fig1)
pp.savefig(fig2)
pp.savefig(fig3)
pp.savefig(fig4)
pp.savefig(fig5)
pp.close()

#Results
pd.ExcelWriter(f'output_{salt_name}_{int(T_filter)}_SI.xlsx', engine='openpyxl') 


#Pitzer
T_filter_df = pd.DataFrame(T_filter, columns=['T Potential'])
Beta_1_MX_MEG_df = pd.DataFrame(Beta_1_MX_MEG, columns=['Beta 1 MX MEG'])
Beta_1_MX_H2O_MEG_df = pd.DataFrame(Beta_1_MX_H2O_MEG, columns=['Beta 1 MX H2O MEG'])

#Chemical Potential
x_MEG_df = pd.DataFrame(x_MEG_SF_All, columns=['x_MEG_SF'])
pot_df = pd.DataFrame(exc_chemical_potential, columns=['Potential exp'])
pot_calc_df = pd.DataFrame(pot_calc, columns=['Potential calc'])
Thetas_pot_df = pd.DataFrame(Thetas_pot, columns=['Thetas pot'])
R_2_pot_df = pd.DataFrame({'R2 pot':[R_2_pot]})
AAD_pot_df = pd.DataFrame({'AAD pot':[AAD_pot]})
maxAD_pot_df = pd.DataFrame({'maxAD pot':[maxAD_pot]})
AARD_pot_df = pd.DataFrame({'AARD pot':[AARD_pot]})
maxARD_pot_df = pd.DataFrame({'maxARD pot':[maxARD_pot]})

#Mixture solubility + Gamma
x_MEG_SF_All_df = pd.DataFrame(x_MEG_SF_All, columns=['x_MEG_sol'])
w_MX_All_df = pd.DataFrame(w_MX_All_True, columns=['w MX'])
b_MX_H2O_MEG_df = pd.DataFrame(b_MX_H2O_MEG, columns=['b_mix exp'])
gamma_MX_H2O_MEG_df = pd.DataFrame(gamma_MX_H2O_MEG, columns=['gamma exp'])
x_MEG_SF_opt_df = pd.DataFrame(x_MEG_SF_opt, columns=['x MEG SF calc'])
w_MEG_SF_opt_df = pd.DataFrame(w_MEG_SF_opt, columns=['w MEG SF calc'])
b_MX_H2O_MEG_optimizated_df = pd.DataFrame(b_MX_H2O_MEG_optimizated, columns=['b_MX_H2O_MEG_optimized'])
gamma_MX_H2O_MEG_optimizated_df = pd.DataFrame(gamma_MX_H2O_MEG_optimizated, columns=['gamma_MX_H2O_MEG_optimized'])
R_2_b_df = pd.DataFrame({'R2 b':[R_2_b]})
AAD_b_df = pd.DataFrame({'AAD b':[AAD_b]})
maxAD_b_df = pd.DataFrame({'maxAD b':[maxAD_b]})
AARD_b_df = pd.DataFrame({'AARD b':[AARD_b]})
maxARD_b_df = pd.DataFrame({'maxARD b':[maxARD_b]})
R_2_g_df = pd.DataFrame({'R2 g':[R_2_gamma]})
AAD_g_df = pd.DataFrame({'AAD g':[AAD_g]})
maxAD_g_df = pd.DataFrame({'maxAD g':[maxAD_g]})
AARD_g_df = pd.DataFrame({'AARD g':[AARD_g]})
maxARD_g_df = pd.DataFrame({'maxARD g':[maxARD_g]})
b_MX_estimado_df = pd.DataFrame(b_MX_H2O_MEG_est, columns=['b_MX_H2O_MEG_est'])
gamma_MX_estimado_df = pd.DataFrame(gamma_MX_H2O_MEG_est, columns=['gamma_MX_H2O_MEG_est'])

result = pd.concat([T_filter_df,x_MEG_df,pot_df, pot_calc_df, Thetas_pot_df, R_2_pot_df, AAD_pot_df, maxAD_pot_df, AARD_pot_df, maxARD_pot_df, Beta_1_MX_MEG_df, Beta_1_MX_H2O_MEG_df, x_MEG_SF_All_df, w_MX_All_df, b_MX_H2O_MEG_df,b_MX_estimado_df, gamma_MX_H2O_MEG_df, gamma_MX_estimado_df, R_2_b_df, AAD_b_df, maxAD_b_df, AARD_b_df, maxARD_b_df, R_2_g_df, AAD_g_df,  maxAD_g_df, AARD_g_df, maxARD_g_df, x_MEG_SF_opt_df, w_MEG_SF_opt_df, b_MX_H2O_MEG_optimizated_df, gamma_MX_H2O_MEG_optimizated_df], axis=1)
result.to_excel(f'output_{salt_name}_{int(T_filter)}_SI.xlsx', sheet_name = salt_name, header=True, index=False)


