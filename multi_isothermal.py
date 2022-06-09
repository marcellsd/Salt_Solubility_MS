# -*- coding: utf-8 -*-
"""
Created on Tue May 25 16:31:14 2021

@author: marcell
"""
import sys

from numpy.core.fromnumeric import mean
sys.path.append('Modules')
from Modules import pitzer
import pandas as pd
from Modules import transformations
from numpy import array, linspace, log, isnan, empty, exp
from Modules import solubility
import matplotlib.pyplot as plt
from Modules import excess
from Modules import statistics
from Modules import dataFilter
from tkinter import filedialog, Tk

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
w_MX_H2O_bin_L, w_H2O_bin_L, T_bin_H2O_L = [], [], []
w_MEG_bin_L, w_MX_MEG_bin_L, T_bin_MEG_L = [], [], []
w_H2O_tern_L, w_MEG_tern_L, w_MX_tern_L, w_H2O_SF_tern_L, w_MEG_SF_tern_L, T_tern_L = [], [], [], [], [], []
w_H2O_All_L, w_MEG_All_L, w_MX_All_L, w_H2O_SF_All_L, w_MEG_SF_All_L, T_All_L = [], [], [], [], [], []   

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
        T_bin_H2O_L.append(T[i])
    elif w_H2O_SF[i] <= 0.01:
        w_MEG_bin_L.append((1-w_H2O_SF[i])*(1 - w_MX_all[i]))
        w_MX_MEG_bin_L.append(w_MX_all[i])
        T_bin_MEG_L.append(T[i])
    else:
        w_H2O_tern_L.append(w_H2O_SF[i]*(1 - w_MX_all[i]))
        w_MX_tern_L.append(w_MX_all[i])
        w_MEG_tern_L.append(1 - w_H2O_SF[i]*(1-w_MX_all[i]) - w_MX_all[i])
        w_H2O_SF_tern_L.append(w_H2O_SF[i])
        w_MEG_SF_tern_L.append(1-w_H2O_SF[i])
        T_tern_L.append(T[i])  
        
#Substitute values of zero
for i in range (0, len(w_MX_MEG_bin_L)):
    if w_MX_MEG_bin_L[i] == 0:
        w_MX_MEG_bin_L.append(1e-6)

#List to array
w_MX_H2O_bin, w_H2O_bin, T_bin_H2O  = array(w_MX_H2O_bin_L), array(w_H2O_bin_L), array(T_bin_H2O_L)
w_MEG_bin, w_MX_MEG_bin, T_bin_MEG  = array(w_MEG_bin_L), array(w_MX_MEG_bin_L), array(T_bin_MEG_L)
w_H2O_tern, w_MEG_tern, w_MX_tern, w_H2O_SF_tern, T_tern = array(w_H2O_tern_L), array(w_MEG_tern_L), array(w_MX_tern_L), array(w_H2O_SF_tern_L), array(T_tern_L)
w_H2O_All, w_MEG_All, w_MX_All, w_H2O_SF_All, w_MEG_SF_All, T_All = array(w_H2O_All_L), array(w_MEG_All_L), array(w_MX_All_L), array(w_H2O_SF_All_L), array(w_MEG_SF_All_L), array(T_All_L)


#Transforming Mass fraction into Molar Salt-Free
x_MEG_SF_All = transformations.x_MEG_SF(w_H2O_SF_All)
x_H2O_SF_All = (1-x_MEG_SF_All)

#Molalities Calculations

b_MX_H2O_bin = transformations.mass_fraction_to_molality(MM_MX, w_MX_H2O_bin)
b_MX_MEG_bin = transformations.mass_fraction_to_molality(MM_MX, w_MX_MEG_bin)

if len(T_filter) != 1:

    Thetas_H2O = solubility.n_linear_regression(T_bin_H2O, b_MX_H2O_bin)
    Thetas_MEG = solubility.n_linear_regression(T_bin_MEG, b_MX_MEG_bin)
    
    T_bin_H2O_a = linspace(min(T_bin_H2O),max(T_bin_H2O), 100, endpoint=True)
    T_bin_MEG_a = linspace(min(T_bin_MEG),max(T_bin_MEG), 100, endpoint=True)
    
    b_H2O_calc = solubility.func_solubility(T_bin_H2O_a, Thetas_H2O)
    b_MEG_calc = solubility.func_solubility(T_bin_MEG_a, Thetas_MEG)
    
    #Calculo do R²

    b_MX_H2O_bin_calc = solubility.func_solubility(T_bin_H2O, Thetas_H2O)
    b_MX_MEG_bin_calc = solubility.func_solubility(T_bin_MEG, Thetas_MEG)
    
    R_2_H2O = statistics.coefficient_of_determination2(b_MX_H2O_bin, b_MX_H2O_bin_calc)
    R_2_MEG = statistics.coefficient_of_determination2(b_MX_MEG_bin, b_MX_MEG_bin_calc)
    
    #Molality x Temperature Graphs
    
    fig0 = plt.figure(0)
    plt.xlabel("T_h2o")
    plt.ylabel("b_h2o")
    plt.plot(T_bin_H2O, b_MX_H2O_bin, 'y.')
    plt.plot(T_bin_H2O_a, b_H2O_calc[:], 'r-')
    plt.title("Model X Experimental h2o")
    
    fig1 = plt.figure(1)
    plt.xlabel("T_MEG")
    plt.ylabel("b_MEG")
    plt.plot(T_bin_MEG, b_MX_MEG_bin, 'y.')
    plt.plot(T_bin_MEG_a, b_MEG_calc[:], 'r-')
    plt.title("Model X Experimental meg")
 

#Mixtures H2O_MEG_MX calculations
if len(T_filter) != 1:
    b_H2O_MX_calc = solubility.func_solubility(T_All, Thetas_H2O)
    b_MEG_MX_calc = solubility.func_solubility(T_All, Thetas_MEG)

else: 
    b_H2O_MX_calc = b_MX_H2O_bin
    b_MEG_MX_calc = b_MX_MEG_bin
  
A_phi_H2O = pitzer.A_phi_H2O(T_All)

rho_H2O, rho_MEG, epsilon_r_H2O, epsilon_r_MEG = pitzer.rho_and_epsilon_r_of_pures(T_All)

A_phi_MEG = pitzer.A_phi_MEG(A_phi_H2O, rho_H2O, epsilon_r_H2O, rho_MEG, epsilon_r_MEG)

rho_H2O_MEG, epsilon_r_H2O_MEG = pitzer.rho_and_epsilon_r_of_mixing(rho_H2O, rho_MEG, epsilon_r_H2O, epsilon_r_MEG, 
                                                                    w_H2O_SF_All, x_MEG_SF_All, x_H2O_SF_All, w_MEG_SF_All, T_All)

A_phi_H2O_MEG = pitzer.A_phi_of_Mixing(rho_H2O_MEG, rho_H2O, epsilon_r_H2O_MEG, epsilon_r_H2O, T_All)

Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O = pitzer.pitzerParameters_MX_H2O(Beta_0, Beta_1, Cphi, T_All, Zm, Rm, Rx)

ln_gamma_MX_H2O = pitzer.ln_gamma_MX(nu_M, nu_X, Zm, Zx, b_H2O_MX_calc, A_phi_H2O, Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O)
Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG = pitzer.pitzerParameter_MX_MEG(nu_M, nu_X, A_phi_MEG, Beta_0_MX_H2O, Beta_1_MX_H2O,
                                                                          Cphi_MX_H2O, ln_gamma_MX_H2O, b_H2O_MX_calc, b_MEG_MX_calc, 
                                                                          T_All, Zm, Zx, Rm, Rx, DeltaG_transf)

Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, Cphi_MX_H2O_MEG = pitzer.Pitzer_parameters_in_mixtures(Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O, 
                                                                                             Beta_1_MX_MEG, epsilon_r_H2O, epsilon_r_MEG, epsilon_r_H2O_MEG)

# Excess chemical potential calculations

b_MX_H2O_MEG = transformations.mass_fraction_to_molality(MM_MX, w_MX_All)
ln_b_MX_H2O_MEG_ideal = excess.ln_b_MX_H2O_MEG_ideal(x_MEG_SF_All ,b_H2O_MX_calc, b_MEG_MX_calc)
ln_gamma_MX_MEG = pitzer.ln_gamma_MX(nu_M, nu_X, Zm, Zx, b_MEG_MX_calc, A_phi_MEG, Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG)
ln_gamma_MX_H2O_MEG_ideal = excess.ln_gamma_MX_H2O_MEG_ideal(x_MEG_SF_All, ln_gamma_MX_H2O, ln_gamma_MX_MEG)
ln_b_MX_H2O_MEG = log(b_MX_H2O_MEG)
ln_b_MX_H2O_MEG_exc = excess.ln_b_MX_H2O_MEG_exc(ln_b_MX_H2O_MEG_ideal, ln_b_MX_H2O_MEG)
ln_gamma_MX_H2O_MEG = pitzer.ln_gamma_MX(nu_M, nu_X, Zm, Zx, b_MX_H2O_MEG, A_phi_H2O_MEG, Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, Cphi_MX_H2O_MEG)  
ln_gamma_MX_H2O_MEG_exc = excess.ln_gamma_MX_H2O_MEG_exc(ln_gamma_MX_H2O_MEG_ideal, ln_gamma_MX_H2O_MEG)
exc_chemical_potential = excess.excess_chemical_potential(ln_b_MX_H2O_MEG_exc, ln_gamma_MX_H2O_MEG_exc)


#Not a Number (NAN) treatment

treated_exc_chem_pot_L = []
T_All_treated_L = []
x_MEG_SF_All_treated_L = []

vector_length_chemical_potential = len(exc_chemical_potential)

for i in range(0,vector_length_chemical_potential):
    if (isnan(exc_chemical_potential[i]) == False):
        treated_exc_chem_pot_L.append(exc_chemical_potential[i])
        T_All_treated_L.append(T_All[i])
        x_MEG_SF_All_treated_L.append(x_MEG_SF_All[i])
  
treated_exc_chem_pot = array(treated_exc_chem_pot_L)
T_All_treated = array(T_All_treated_L)
x_MEG_SF_All_treated = array(x_MEG_SF_All_treated_L)

#Number of parameters in chemical Potential model (choose 4 or 6)
n = 6
     
Thetas_pot = excess.n_linear_regression_chemical_potential(treated_exc_chem_pot, T_All_treated, x_MEG_SF_All_treated, n)

pot_calc = excess.model_excess_chemical_potential(T_All_treated, x_MEG_SF_All_treated, Thetas_pot, n)

x_MEG = linspace( 0.0,1.0, 100, endpoint=True)

t_filter_length = len(T_filter)
pot_calc_matrix = empty([100, t_filter_length])

for i in range(0,t_filter_length):
    pot_calc_matrix[:,i] = excess.model_excess_chemical_potential(T_filter[i], x_MEG, Thetas_pot, n)

#Create lists
pot_calc_L = dataFilter.filter_1(treated_exc_chem_pot, T_All_treated,T_filter)
x_MEG_SF_All_treated_L =dataFilter.filter_1(x_MEG_SF_All_treated,T_All_treated,T_filter)

#Chemical Potential statistics
R_2_pot = statistics.coefficient_of_determination2(treated_exc_chem_pot, pot_calc)
AAD_pot, maxAD_pot, AARD_pot, maxARD_pot = statistics.DAM_DMR(treated_exc_chem_pot, pot_calc)

#Filtering pitzer parameters by temperature
Beta_0_MX_H2O_filter, Beta_1_MX_H2O_filter, Cphi_MX_H2O_filter = dataFilter.filter(Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O, T_All, T_filter)
Beta_0_MX_MEG_filter, Beta_1_MX_MEG_filter, Cphi_MX_MEG_filter = dataFilter.filter(Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG, T_All, T_filter)

#Mixture solubility calculations
number_of_points = 31

x_MEG_SF_optimizated = linspace(0.0,1.0,number_of_points,endpoint=True)
x_H2O_SF_optimizated = 1-x_MEG_SF_optimizated
w_MEG_SF_optimizated = transformations.w_MEG_SF(x_MEG_SF_optimizated)
w_H2O_SF_optimizated = 1-w_MEG_SF_optimizated

b_MX_H2O_MEG_optimizated = empty([number_of_points, t_filter_length])
gamma_MX_H2O_MEG_optimizated = empty([number_of_points, t_filter_length])
for i in range(0, t_filter_length):
    b_MX_H2O_MEG_optimizated[:,i], gamma_MX_H2O_MEG_optimizated[:,i] = excess.optimization_b_MX_H2O_MEG(x_MEG_SF_optimizated, w_H2O_SF_optimizated, x_H2O_SF_optimizated, w_MEG_SF_optimizated, T_filter[i], 
                                                                                                     Thetas_pot, Thetas_H2O, Thetas_MEG, 
                                                                                                     Beta_0_MX_H2O_filter[i], Beta_1_MX_H2O_filter[i], Cphi_MX_H2O_filter[i],
                                                                                                     Beta_0_MX_MEG_filter[i], Beta_1_MX_MEG_filter[i], Cphi_MX_MEG_filter[i],
                                                                                                     Zm, nu_M, nu_X, Zx, n)


x_MEG_SF_All_L =dataFilter.filter_1(x_MEG_SF_All,T_All,T_filter)
b_MX_H2O_MEG_L = dataFilter.filter_1(b_MX_H2O_MEG, T_All,T_filter)
gamma_MX_H2O_MEG = exp(ln_gamma_MX_H2O_MEG)
gamma_MX_H2O_MEG_L = dataFilter.filter_1(gamma_MX_H2O_MEG, T_All,T_filter)

#Prediction calculations
T_pred = array([330.0, 400.0])
pot_calc_pred = empty([100, len(T_pred)])
for i in range(0,len(T_pred)):
    pot_calc_pred[:,i] = excess.model_excess_chemical_potential(T_pred[i], x_MEG, Thetas_pot, n)
#Pitzer predictions calculations
A_phi_H2O_pred = pitzer.A_phi_H2O(T_pred)
rho_H2O_pred, rho_MEG_pred, epsilon_r_H2O_pred, epsilon_r_MEG_pred = pitzer.rho_and_epsilon_r_of_pures(T_pred)
A_phi_MEG_pred = pitzer.A_phi_MEG(A_phi_H2O_pred, rho_H2O_pred, epsilon_r_H2O_pred, rho_MEG_pred, epsilon_r_MEG_pred)
Beta_0_MX_H2O_pred, Beta_1_MX_H2O_pred, Cphi_MX_H2O_pred = pitzer.pitzerParameters_MX_H2O(Beta_0, Beta_1, Cphi, T_pred, Zm, Rm, Rx)

b_H2O_MX_calc_pred = solubility.func_solubility(T_pred, Thetas_H2O)
b_MEG_MX_calc_pred = solubility.func_solubility(T_pred, Thetas_MEG)

ln_gamma_MX_H2O_pred = pitzer.ln_gamma_MX(nu_M, nu_X, Zm, Zx, b_H2O_MX_calc_pred, A_phi_H2O_pred, Beta_0_MX_H2O_pred, Beta_1_MX_H2O_pred, Cphi_MX_H2O_pred)
Beta_0_MX_MEG_pred, Beta_1_MX_MEG_pred, Cphi_MX_MEG_pred = pitzer.pitzerParameter_MX_MEG(nu_M, nu_X, A_phi_MEG_pred, Beta_0_MX_H2O_pred, Beta_1_MX_H2O_pred,
                                                                          Cphi_MX_H2O_pred, ln_gamma_MX_H2O_pred, b_H2O_MX_calc_pred, b_MEG_MX_calc_pred, 
                                                                          T_pred, Zm, Zx, Rm, Rx, DeltaG_transf)


b_MX_H2O_MEG_pred = empty([number_of_points, len(T_pred)])
gamma_MX_H2O_MEG_pred = empty([number_of_points, len(T_pred)])
for i in range(0, len(T_pred)):
    b_MX_H2O_MEG_pred[:,i], gamma_MX_H2O_MEG_pred[:,i] = excess.optimization_b_MX_H2O_MEG(x_MEG_SF_optimizated, w_H2O_SF_optimizated, x_H2O_SF_optimizated, w_MEG_SF_optimizated, T_pred[i], 
                                                                                                     Thetas_pot, Thetas_H2O, Thetas_MEG, 
                                                                                                     Beta_0_MX_H2O_pred[i], Beta_1_MX_H2O_pred[i], Cphi_MX_H2O_pred[i],
                                                                                                     Beta_0_MX_MEG_pred[i], Beta_1_MX_MEG_pred[i], Cphi_MX_MEG_pred[i],
                                                                                                     Zm, nu_M, nu_X, Zx, n)

b_MX_H2O_MEG_est = empty([number_of_points, len(T_All)])
gamma_MX_H2O_MEG_est = empty([number_of_points, len(T_All)])
b_MX_H2O_MEG_est, gamma_MX_H2O_MEG_est = excess.optimization_b_MX_H2O_MEG_statistics(x_MEG_SF_All, w_H2O_SF_All, x_H2O_SF_All, w_MEG_SF_All, T_All, 
                                                                                                     Thetas_pot, Thetas_H2O, Thetas_MEG, 
                                                                                                     Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O,
                                                                                                     Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG,
                                                                                                     Zm, nu_M, nu_X, Zx, n)
#Solubility and Gamma Statistics
R_2_b=statistics.coefficient_of_determination2(b_MX_H2O_MEG, b_MX_H2O_MEG_est)
R_2_gamma=statistics.coefficient_of_determination2(gamma_MX_H2O_MEG, gamma_MX_H2O_MEG_est)
AAD_b, maxAD_b, AARD_b, maxARD_b = statistics.DAM_DMR(b_MX_H2O_MEG, b_MX_H2O_MEG_est)
AAD_g, maxAD_g, AARD_g, maxARD_g = statistics.DAM_DMR(gamma_MX_H2O_MEG, gamma_MX_H2O_MEG_est)

#Graphics 
#Chemical Potential Graphics
number_of_plots = len(T_filter)
colormap = plt.cm.tab20b
colors = [colormap(i) for i in linspace(0, 1,number_of_plots)] 

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
ax2.set_prop_cycle('color', colors)        
plt.figure(2)
plt.xlabel("x_MEG_SF")
plt.ylabel("Potential")
for i in range(0,t_filter_length):
    plt.plot(x_MEG_SF_All_treated_L[i], pot_calc_L[i], '.')
for i in range(0, t_filter_length):
    plt.plot(x_MEG, pot_calc_matrix[:,i], label=T_filter[i])
plt.title("Calculated Chemical Potentials")
plt.legend()

#Graphic b
fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111)
ax3.set_prop_cycle('color', colors)
plt.xlabel("x_MEG_SF")
plt.ylabel("b")
for i in range(0,t_filter_length):
    plt.plot(x_MEG_SF_All_L[i], b_MX_H2O_MEG_L[i], '.')
for i in range(0, t_filter_length):
    plt.plot(x_MEG_SF_optimizated, b_MX_H2O_MEG_optimizated[:,i], label=T_filter[i])
plt.title("Calculated b_MX_H2O_MEG")
plt.legend()

#Graphic gamma
fig4 = plt.figure(4)
ax4 = fig4.add_subplot(111)
ax4.set_prop_cycle('color', colors)
plt.xlabel("x_MEG_SF")
plt.ylabel("gamma")
for i in range(0,t_filter_length):
    plt.plot(x_MEG_SF_All_L[i],gamma_MX_H2O_MEG_L[i], '.')
for i in range(0, t_filter_length):
    plt.plot(x_MEG_SF_optimizated, gamma_MX_H2O_MEG_optimizated[:,i], label=T_filter[i])
plt.title("Calculated gamma_MX_H2O_MEG ")
plt.legend()

fig5 = plt.figure(5)
ax5 = fig5.add_subplot(111)
ax5.set_prop_cycle('color', colors)
plt.xlabel("x_MEG_SF")
plt.ylabel("Potential")
for i in range(0,len(T_pred)):
    plt.plot(x_MEG, pot_calc_pred[:,i], label=T_pred[i])
plt.title("Predicted Chemical Potential")
plt.legend()

fig6 = plt.figure(6)
ax6 = fig6.add_subplot(111)
ax6.set_prop_cycle('color', colors)
plt.xlabel("x_MEG_SF")
plt.ylabel("b")
for i in range(0, len(T_pred)):
    plt.plot(x_MEG_SF_optimizated, b_MX_H2O_MEG_pred[:,i], label=T_pred[i])
plt.title("Predicted b_MX_H2O_MEG")
plt.legend()

fig7 = plt.figure(7)
ax7 = fig7.add_subplot(111)
ax7.set_prop_cycle('color', colors)
plt.xlabel("x_MEG_SF")
plt.ylabel("gamma")
for i in range(0, len(T_pred)):
    plt.plot(x_MEG_SF_optimizated, gamma_MX_H2O_MEG_pred[:,i], label=T_pred[i])
plt.title("Predicted gamma_MX_H2O_MEG")
plt.legend()


from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages(f'output_{salt_name}_MI.pdf')
plt.savefig(pp, format='pdf')
pp.savefig(fig0)
pp.savefig(fig1)
pp.savefig(fig2)
pp.savefig(fig3)
pp.savefig(fig4)
pp.savefig(fig5)
pp.savefig(fig6)
pp.savefig(fig7)
pp.close()

#Results
pd.ExcelWriter(f'output_{salt_name}_MI.xlsx', engine='openpyxl') 


#Pure components solubility
T_bin_H2O_df = pd.DataFrame(T_bin_H2O, columns=['T_H2O'])
T_bin_MEG_df = pd.DataFrame(T_bin_MEG, columns=['T_MEG'])
b_bin_H2O_df = pd.DataFrame(b_MX_H2O_bin, columns=['b_MX_H2O'])
b_bin_MEG_df = pd.DataFrame(b_MX_MEG_bin, columns=['b_MX_MEG']) 
Thetas_H2O_df = pd.DataFrame(Thetas_H2O, columns=['Thetas H2O'])
Thetas_MEG_df = pd.DataFrame(Thetas_MEG, columns=['Thetas MEG'])
r2_H2O_df = pd.DataFrame({'R2 H2O':[R_2_H2O]})
r2_MEG_df = pd.DataFrame({'R2 MEG':[R_2_MEG]})

#Pitzer
Beta_1_MX_MEG_df = pd.DataFrame(Beta_1_MX_MEG, columns=['Beta 1 MX MEG'])
Beta_1_MX_H2O_MEG_df = pd.DataFrame(Beta_1_MX_H2O_MEG, columns=['Beta 1 MX H2O MEG'])
T_all_df1 = pd.DataFrame(T_All, columns=['T all'])

#Chemical potential
x_MEG_df = pd.DataFrame(x_MEG_SF_All_treated, columns=['x_MEG_SF'])
pot_df = pd.DataFrame(treated_exc_chem_pot, columns=['Potential exp'])
pot_calc_df = pd.DataFrame(pot_calc, columns=['Potential calc'])
T_All_treated_df = pd.DataFrame(T_All_treated, columns=['T Potential'])
Thetas_pot_df = pd.DataFrame(Thetas_pot, columns=['Thetas pot'])
R_2_pot_df = pd.DataFrame({'R2 pot':[R_2_pot]})
AAD_pot_df = pd.DataFrame({'AAD pot':[AAD_pot]})
maxAD_pot_df = pd.DataFrame({'maxAD pot':[maxAD_pot]})
AARD_pot_df = pd.DataFrame({'AARD pot':[AARD_pot]})
maxARD_pot_df = pd.DataFrame({'maxARD pot':[maxARD_pot]})

#Mixture solubility + Gamma
x_MEG_SF_All_df = pd.DataFrame(x_MEG_SF_All, columns=['x_MEG_sol'])
w_MX_All_df = pd.DataFrame(w_MX_All, columns=['w MX'])
T_all_df2 = pd.DataFrame(T_All, columns=['T sol mix'])
b_MX_H2O_MEG_df = pd.DataFrame(b_MX_H2O_MEG, columns=['b_mix exp'])
gamma_MX_H2O_MEG_df = pd.DataFrame(gamma_MX_H2O_MEG, columns=['gamma exp'])
x_MEG_SF_optimizated_df = pd.DataFrame(x_MEG_SF_optimizated, columns=['x MEG SF calc'])
w_MEG_SF_optimizated_df = pd.DataFrame(w_MEG_SF_optimizated, columns=['w MEG SF calc'])
b_MX_H2O_MEG_optimizated_df = pd.DataFrame(b_MX_H2O_MEG_optimizated)
gamma_MX_H2O_MEG_optimizated_df = pd.DataFrame(gamma_MX_H2O_MEG_optimizated)
T_filter_df = pd.DataFrame(T_filter, columns=['Isothermal Temperatures'])
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
b_MX_estimated_df = pd.DataFrame(b_MX_H2O_MEG_est, columns=['b_MX_H2O_MEG_est'])
gamma_MX_estimated_df = pd.DataFrame(gamma_MX_H2O_MEG_est, columns=['gamma_MX_H2O_MEG_est'])


result = pd.concat([b_bin_H2O_df,T_bin_H2O_df, b_bin_MEG_df, T_bin_MEG_df,Thetas_H2O_df,Thetas_MEG_df,r2_H2O_df,r2_MEG_df,x_MEG_df,pot_df, pot_calc_df, T_All_treated_df, Thetas_pot_df, R_2_pot_df, AAD_pot_df, maxAD_pot_df, AARD_pot_df, maxARD_pot_df, Beta_1_MX_MEG_df, Beta_1_MX_H2O_MEG_df, T_all_df1, x_MEG_SF_All_df, w_MX_All_df,T_all_df2, b_MX_H2O_MEG_df,b_MX_estimated_df, gamma_MX_H2O_MEG_df, gamma_MX_estimated_df, R_2_b_df, AAD_b_df, maxAD_b_df, AARD_b_df, maxARD_b_df, R_2_g_df, AAD_g_df,  maxAD_g_df, AARD_g_df, maxARD_g_df, x_MEG_SF_optimizated_df, w_MEG_SF_optimizated_df, T_filter_df, b_MX_H2O_MEG_optimizated_df, gamma_MX_H2O_MEG_optimizated_df], axis=1)
result.to_excel(f'output_{salt_name}_MI.xlsx', sheet_name = salt_name, header=True, index=False)

