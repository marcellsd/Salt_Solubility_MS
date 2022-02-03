# -*- coding: utf-8 -*-
#%%
"""
Created on Tue May 25 16:31:14 2021

@author: marcell
"""
import sys

from numpy.core.fromnumeric import mean
sys.path.append('Modulos')
import platform
from Modulos import pitzer
#from Modulos import seleciona
import pandas as pd
from Modulos import transformacoes
from numpy import array, linspace, log, isnan, empty, exp
from Modulos import solubilidade
import matplotlib.pyplot as plt
from Modulos import excesso
from Modulos import estatistica
from Modulos import filtro_dados

from tkinter import filedialog, Tk

root = Tk()
root.withdraw()
data_file_name = filedialog.askopenfilename()

df = pd.read_excel(data_file_name, index_col=('Parametro'))
Rm = df.loc['Raio cation', 'ParamValor']
Rx = df.loc['Raio anion', 'ParamValor']
Zm = df.loc['Carga cation', 'ParamValor']
Zx = df.loc['Carga anion', 'ParamValor']
nu_M = df.loc['Coef esteq cation', 'ParamValor']
nu_X = df.loc['Coef esteq anion', 'ParamValor']
mm_MX = df.loc['Massa molar sal (g/mol)', 'ParamValor']
MM_MX = mm_MX/1000 #Kg/mol
DeltaG_transf = df.loc['DeltaG transf (J/mol)', 'ParamValor']
nome_cation = df.loc['Cation', 'ParamValor']
nome_anion = df.loc['Anion', 'ParamValor']

if nu_M==1 and nu_X==1:
    nome_sal = nome_cation + nome_anion
elif nu_M==1 and nu_X!=1:
    nome_sal = nome_cation + nome_anion + str(nu_X)
elif nu_M!=1 and nu_X==1:
    nome_sal = nome_cation + str(nu_M) + nome_anion
else:
    nome_sal = nome_cation + str(nu_M) + nome_anion + str(nu_X)

#Coef. Pitzer
Beta_0 = df.loc['Beta0_ref_25C', 'ParamValor']
Beta_1 = df.loc['Beta1_ref_25C', 'ParamValor']
Cphi = df.loc['Cphi_ref_25C', 'ParamValor']

data_file_sol = pd.read_excel(data_file_name, 'dados')
T = pd.DataFrame(data_file_sol, columns=['T (K)']).to_numpy()
w_H2O_SF = pd.DataFrame(data_file_sol, columns=['w_H20_SF']).to_numpy()
w_MX_Todos = pd.DataFrame(data_file_sol, columns=['wsal']).to_numpy()
T_filtro = pd.DataFrame(data_file_sol, columns=['T (K)']).drop_duplicates().to_numpy()

#Dados que w_H2O_SF ≃ 1 -> Binários Água-Sal
#Dados que w_H2O_SF ≃ 0 -> Binários MEG-Sal

#Separação dos dados em ternário e binários
w_MX_H2O_bin_L, w_H2O_bin_L, T_bin_H2O_L = [], [], []

w_MEG_bin_L, w_MX_MEG_bin_L, T_bin_MEG_L = [], [], []

w_H2O_tern_L, w_MEG_tern_L, w_MX_tern_L, w_H2O_SF_tern_L, w_MEG_SF_tern_L, T_tern_L = [], [], [], [], [], []

w_H2O_All_L, w_MEG_All_L, w_MX_All_L, w_H2O_SF_All_L, w_MEG_SF_All_L, T_All_L = [], [], [], [], [], []   

tamanho_vetor = len(w_H2O_SF)
for i in range (0, tamanho_vetor):
    w_H2O_All_L.append(w_H2O_SF[i]*(1 - w_MX_Todos[i]))
    w_MX_All_L.append(w_MX_Todos[i])
    w_MEG_All_L.append(1 - w_H2O_SF[i]*(1-w_MX_Todos[i]) - w_MX_Todos[i])
    w_H2O_SF_All_L.append(w_H2O_SF[i])
    w_MEG_SF_All_L.append(1-w_H2O_SF[i])
    T_All_L.append(T[i])
    
    if w_H2O_SF[i] >= 0.99:
        w_MX_H2O_bin_L.append(w_MX_Todos[i])
        w_H2O_bin_L.append(w_H2O_SF[i]*(1 - w_MX_Todos[i]))
        T_bin_H2O_L.append(T[i])
    elif w_H2O_SF[i] <= 0.01:
        w_MEG_bin_L.append((1-w_H2O_SF[i])*(1 - w_MX_Todos[i]))
        w_MX_MEG_bin_L.append(w_MX_Todos[i])
        T_bin_MEG_L.append(T[i])
    else:
        w_H2O_tern_L.append(w_H2O_SF[i]*(1 - w_MX_Todos[i]))
        w_MX_tern_L.append(w_MX_Todos[i])
        w_MEG_tern_L.append(1 - w_H2O_SF[i]*(1-w_MX_Todos[i]) - w_MX_Todos[i])
        w_H2O_SF_tern_L.append(w_H2O_SF[i])
        w_MEG_SF_tern_L.append(1-w_H2O_SF[i])
        T_tern_L.append(T[i])  
        
#Substituir valores 0
for i in range (0, len(w_MX_MEG_bin_L)):
    if w_MX_MEG_bin_L[i] == 0:
        w_MX_MEG_bin_L.append(1e-6)

#Transformando de lista para array
w_MX_H2O_bin, w_H2O_bin, T_bin_H2O  = array(w_MX_H2O_bin_L), array(w_H2O_bin_L), array(T_bin_H2O_L)

w_MEG_bin, w_MX_MEG_bin, T_bin_MEG  = array(w_MEG_bin_L), array(w_MX_MEG_bin_L), array(T_bin_MEG_L)

w_H2O_tern, w_MEG_tern, w_MX_tern, w_H2O_SF_tern, T_tern = array(w_H2O_tern_L), array(w_MEG_tern_L), array(w_MX_tern_L), array(w_H2O_SF_tern_L), array(T_tern_L)

w_H2O_All, w_MEG_All, w_MX_All, w_H2O_SF_All, w_MEG_SF_All, T_All = array(w_H2O_All_L), array(w_MEG_All_L), array(w_MX_All_L), array(w_H2O_SF_All_L), array(w_MEG_SF_All_L), array(T_All_L)

#Transformando Fração Mássica para Molar SF

x_MEG_SF_All = transformacoes.x_MEG_SF(w_H2O_SF_All)
x_H2O_SF_All = (1-x_MEG_SF_All)

#Calculo das molalidades

b_MX_H2O_bin = transformacoes.fracao_massica_para_molalidade(MM_MX, w_MX_H2O_bin)
b_MX_MEG_bin = transformacoes.fracao_massica_para_molalidade(MM_MX, w_MX_MEG_bin)

if len(T_filtro) != 1:

    Thetas_H2O = solubilidade.regressao_n_linear(T_bin_H2O, b_MX_H2O_bin)
    Thetas_MEG = solubilidade.regressao_n_linear(T_bin_MEG, b_MX_MEG_bin)
    
    
    T_bin_H2O_a = linspace(min(T_bin_H2O),max(T_bin_H2O), 100, endpoint=True)
    T_bin_MEG_a = linspace(min(T_bin_MEG),max(T_bin_MEG), 100, endpoint=True)
    
    b_H2O_calc = solubilidade.func_solubilidade(T_bin_H2O_a, Thetas_H2O)
    b_MEG_calc = solubilidade.func_solubilidade(T_bin_MEG_a, Thetas_MEG)
    
    #Calculo do R²
    b_MX_H2O_bin_calc = solubilidade.func_solubilidade(T_bin_H2O, Thetas_H2O)
    b_MX_MEG_bin_calc = solubilidade.func_solubilidade(T_bin_MEG, Thetas_MEG)
    
    R_2_H2O = estatistica.coef_determinacao2(b_MX_H2O_bin, b_MX_H2O_bin_calc)
    R_2_MEG = estatistica.coef_determinacao2(b_MX_MEG_bin, b_MX_MEG_bin_calc)
    
    #Gráficos Molalidade x Temperatura
    
    fig0 = plt.figure(0)
    plt.xlabel("T_h2o")
    plt.ylabel("b_h2o")
    plt.plot(T_bin_H2O, b_MX_H2O_bin, 'y.')
    plt.plot(T_bin_H2O_a, b_H2O_calc[:], 'r-')
    plt.title("Modelo X Experimental h2o")
    
    
    fig1 = plt.figure(1)
    plt.xlabel("T_MEG")
    plt.ylabel("b_MEG")
    plt.plot(T_bin_MEG, b_MX_MEG_bin, 'y.')
    plt.plot(T_bin_MEG_a, b_MEG_calc[:], 'r-')
    plt.title("Modelo X Experimental meg")
 

#Cálculo para a mistura H2O_MEG_MX
if len(T_filtro) != 1:
    b_H2O_MX_calc = solubilidade.func_solubilidade(T_All, Thetas_H2O)
    b_MEG_MX_calc = solubilidade.func_solubilidade(T_All, Thetas_MEG)
else: 
    b_H2O_MX_calc = b_MX_H2O_bin
    b_MEG_MX_calc = b_MX_MEG_bin
    #print(f'bH2O: {b_H2O_MX_calc}')
    #print(f'bMEG: {b_MEG_MX_calc}')
  
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

#Cálculo do Potencial Químico de excesso

b_MX_H2O_MEG = transformacoes.fracao_massica_para_molalidade(MM_MX, w_MX_All)
ln_b_MX_H2O_MEG_ideal = excesso.ln_b_MX_H2O_MEG_ideal(x_MEG_SF_All ,b_H2O_MX_calc, b_MEG_MX_calc)
ln_gamma_MX_MEG = pitzer.ln_gamma_MX(nu_M, nu_X, Zm, Zx, b_MEG_MX_calc, A_phi_MEG, Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG)
ln_gamma_MX_H2O_MEG_ideal = excesso.ln_gamma_MX_H2O_MEG_ideal(x_MEG_SF_All, ln_gamma_MX_H2O, ln_gamma_MX_MEG)
ln_b_MX_H2O_MEG = log(b_MX_H2O_MEG)
ln_b_MX_H2O_MEG_exc = excesso.ln_b_MX_H2O_MEG_exc(ln_b_MX_H2O_MEG_ideal, ln_b_MX_H2O_MEG)
ln_gamma_MX_H2O_MEG = pitzer.ln_gamma_MX(nu_M, nu_X, Zm, Zx, b_MX_H2O_MEG, A_phi_H2O_MEG, Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, Cphi_MX_H2O_MEG)  
ln_gamma_MX_H2O_MEG_exc = excesso.ln_gamma_MX_H2O_MEG_exc(ln_gamma_MX_H2O_MEG_ideal, ln_gamma_MX_H2O_MEG)
potencial_quimico_exc = excesso.potencial_quimico_exc(ln_b_MX_H2O_MEG_exc, ln_gamma_MX_H2O_MEG_exc)
#%%
#Tratar NAN

potencial_quimico_exc_tratado_L = []
T_All_tratado_L = []
x_MEG_SF_All_tratado_L = []

tamanho_vetor_potencial_quimico = len(potencial_quimico_exc)

for i in range(0,tamanho_vetor_potencial_quimico):
    if (isnan(potencial_quimico_exc[i]) == False):
        potencial_quimico_exc_tratado_L.append(potencial_quimico_exc[i])
        T_All_tratado_L.append(T_All[i])
        x_MEG_SF_All_tratado_L.append(x_MEG_SF_All[i])
  
potencial_quimico_exc_tratado = array(potencial_quimico_exc_tratado_L)
T_All_tratado = array(T_All_tratado_L)
x_MEG_SF_All_tratado = array(x_MEG_SF_All_tratado_L)

#Quantidade de parametros potencial quimico (escolher 4 ou 6)
n = 6
     
Thetas_pot = excesso.regressao_n_linear_potencial_quimico(potencial_quimico_exc_tratado, T_All_tratado, x_MEG_SF_All_tratado, n)

pot_calc = excesso.func_potencial_quimico_exc(T_All_tratado, x_MEG_SF_All_tratado, Thetas_pot, n)

x_MEG = linspace( 0.0,1.0, 100, endpoint=True)

tamanho_t_filtro = len(T_filtro)
pot_calc_matrix = empty([100, tamanho_t_filtro])

for i in range(0,tamanho_t_filtro):
    pot_calc_matrix[:,i] = excesso.func_potencial_quimico_exc(T_filtro[i], x_MEG, Thetas_pot, n)

#Criar listas
pot_calc_L = filtro_dados.filtro_1(potencial_quimico_exc_tratado, T_All_tratado,T_filtro)
x_MEG_SF_All_tratado_L =filtro_dados.filtro_1(x_MEG_SF_All_tratado,T_All_tratado,T_filtro)

#Estatistica potencial
R_2_pot = estatistica.coef_determinacao2(potencial_quimico_exc_tratado, pot_calc)
AAD_pot, maxAD_pot, AARD_pot, maxARD_pot = estatistica.DAM_DMR(potencial_quimico_exc_tratado, pot_calc)

#Filtro parametros pitzer por temperatura
Beta_0_MX_H2O_filtro, Beta_1_MX_H2O_filtro, Cphi_MX_H2O_filtro = filtro_dados.filtro(Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O, T_All, T_filtro)
Beta_0_MX_MEG_filtro, Beta_1_MX_MEG_filtro, Cphi_MX_MEG_filtro = filtro_dados.filtro(Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG, T_All, T_filtro)

#Calculo solubilidade Mistura
numero_de_pontos = 31

x_MEG_SF_oti = linspace(0.0,1.0,numero_de_pontos,endpoint=True)
x_H2O_SF_oti = 1-x_MEG_SF_oti
w_MEG_SF_oti = transformacoes.w_MEG_SF(x_MEG_SF_oti)
w_H2O_SF_oti = 1-w_MEG_SF_oti

b_MX_H2O_MEG_otimizacao = empty([numero_de_pontos, tamanho_t_filtro])
gamma_MX_H2O_MEG_otimizacao = empty([numero_de_pontos, tamanho_t_filtro])
for i in range(0, tamanho_t_filtro):
    b_MX_H2O_MEG_otimizacao[:,i], gamma_MX_H2O_MEG_otimizacao[:,i] = excesso.otimizacao_b_MX_H2O_MEG(x_MEG_SF_oti, w_H2O_SF_oti, x_H2O_SF_oti, w_MEG_SF_oti, T_filtro[i], 
                                                                                                     Thetas_pot, Thetas_H2O, Thetas_MEG, 
                                                                                                     Beta_0_MX_H2O_filtro[i], Beta_1_MX_H2O_filtro[i], Cphi_MX_H2O_filtro[i],
                                                                                                     Beta_0_MX_MEG_filtro[i], Beta_1_MX_MEG_filtro[i], Cphi_MX_MEG_filtro[i],
                                                                                                     Zm, nu_M, nu_X, Zx, n)


x_MEG_SF_All_L =filtro_dados.filtro_1(x_MEG_SF_All,T_All,T_filtro)
b_MX_H2O_MEG_L = filtro_dados.filtro_1(b_MX_H2O_MEG, T_All,T_filtro)
gamma_MX_H2O_MEG = exp(ln_gamma_MX_H2O_MEG)
gamma_MX_H2O_MEG_L = filtro_dados.filtro_1(gamma_MX_H2O_MEG, T_All,T_filtro)

#Calculos Predição
T_pred = array([330.0, 400.0])
pot_calc_pred = empty([100, len(T_pred)])
for i in range(0,len(T_pred)):
    pot_calc_pred[:,i] = excesso.func_potencial_quimico_exc(T_pred[i], x_MEG, Thetas_pot, n)
#Pitzer pred
A_phi_H2O_pred = pitzer.A_phi_H2O(T_pred)
rho_H2O_pred, rho_MEG_pred, epsilon_r_H2O_pred, epsilon_r_MEG_pred = pitzer.rho_and_epsilon_r_of_pures(T_pred)
A_phi_MEG_pred = pitzer.A_phi_MEG(A_phi_H2O_pred, rho_H2O_pred, epsilon_r_H2O_pred, rho_MEG_pred, epsilon_r_MEG_pred)
Beta_0_MX_H2O_pred, Beta_1_MX_H2O_pred, Cphi_MX_H2O_pred = pitzer.pitzerParameters_MX_H2O(Beta_0, Beta_1, Cphi, T_pred, Zm, Rm, Rx)

b_H2O_MX_calc_pred = solubilidade.func_solubilidade(T_pred, Thetas_H2O)
b_MEG_MX_calc_pred = solubilidade.func_solubilidade(T_pred, Thetas_MEG)

ln_gamma_MX_H2O_pred = pitzer.ln_gamma_MX(nu_M, nu_X, Zm, Zx, b_H2O_MX_calc_pred, A_phi_H2O_pred, Beta_0_MX_H2O_pred, Beta_1_MX_H2O_pred, Cphi_MX_H2O_pred)
Beta_0_MX_MEG_pred, Beta_1_MX_MEG_pred, Cphi_MX_MEG_pred = pitzer.pitzerParameter_MX_MEG(nu_M, nu_X, A_phi_MEG_pred, Beta_0_MX_H2O_pred, Beta_1_MX_H2O_pred,
                                                                          Cphi_MX_H2O_pred, ln_gamma_MX_H2O_pred, b_H2O_MX_calc_pred, b_MEG_MX_calc_pred, 
                                                                          T_pred, Zm, Zx, Rm, Rx, DeltaG_transf)


b_MX_H2O_MEG_pred = empty([numero_de_pontos, len(T_pred)])
gamma_MX_H2O_MEG_pred = empty([numero_de_pontos, len(T_pred)])
for i in range(0, len(T_pred)):
    b_MX_H2O_MEG_pred[:,i], gamma_MX_H2O_MEG_pred[:,i] = excesso.otimizacao_b_MX_H2O_MEG(x_MEG_SF_oti, w_H2O_SF_oti, x_H2O_SF_oti, w_MEG_SF_oti, T_pred[i], 
                                                                                                     Thetas_pot, Thetas_H2O, Thetas_MEG, 
                                                                                                     Beta_0_MX_H2O_pred[i], Beta_1_MX_H2O_pred[i], Cphi_MX_H2O_pred[i],
                                                                                                     Beta_0_MX_MEG_pred[i], Beta_1_MX_MEG_pred[i], Cphi_MX_MEG_pred[i],
                                                                                                     Zm, nu_M, nu_X, Zx, n)

b_MX_H2O_MEG_est = empty([numero_de_pontos, len(T_All)])
gamma_MX_H2O_MEG_est = empty([numero_de_pontos, len(T_All)])
b_MX_H2O_MEG_est, gamma_MX_H2O_MEG_est = excesso.otimizacao_b_MX_H2O_MEG_estatistica(x_MEG_SF_All, w_H2O_SF_All, x_H2O_SF_All, w_MEG_SF_All, T_All, 
                                                                                                     Thetas_pot, Thetas_H2O, Thetas_MEG, 
                                                                                                     Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O,
                                                                                                     Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG,
                                                                                                     Zm, nu_M, nu_X, Zx, n)
#Estatistica solubilidade e gamma
R_2_b=estatistica.coef_determinacao2(b_MX_H2O_MEG, b_MX_H2O_MEG_est)
R_2_gamma=estatistica.coef_determinacao2(gamma_MX_H2O_MEG, gamma_MX_H2O_MEG_est)
AAD_b, maxAD_b, AARD_b, maxARD_b = estatistica.DAM_DMR(b_MX_H2O_MEG, b_MX_H2O_MEG_est)
AAD_g, maxAD_g, AARD_g, maxARD_g = estatistica.DAM_DMR(gamma_MX_H2O_MEG, gamma_MX_H2O_MEG_est)

#Graficos 
#Grafico Potencial quimico
number_of_plots = len(T_filtro)
colormap = plt.cm.tab20b
colors = [colormap(i) for i in linspace(0, 1,number_of_plots)] 

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)
ax2.set_prop_cycle('color', colors)        
plt.figure(2)
plt.xlabel("x_MEG_SF")
plt.ylabel("Potencial")
for i in range(0,tamanho_t_filtro):
    plt.plot(x_MEG_SF_All_tratado_L[i], pot_calc_L[i], '.')
for i in range(0, tamanho_t_filtro):
    plt.plot(x_MEG, pot_calc_matrix[:,i], label=T_filtro[i])
plt.title("Potenciais Quimicos Calc")
plt.legend()

#Grafico b
fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111)
ax3.set_prop_cycle('color', colors)
plt.xlabel("x_MEG_SF")
plt.ylabel("b")
for i in range(0,tamanho_t_filtro):
    plt.plot(x_MEG_SF_All_L[i], b_MX_H2O_MEG_L[i], '.')
for i in range(0, tamanho_t_filtro):
    plt.plot(x_MEG_SF_oti, b_MX_H2O_MEG_otimizacao[:,i], label=T_filtro[i])
plt.title("b_MX_H2O_MEG Calc")
plt.legend()

#Grafico gamma
fig4 = plt.figure(4)
ax4 = fig4.add_subplot(111)
ax4.set_prop_cycle('color', colors)
plt.xlabel("x_MEG_SF")
plt.ylabel("gamma")
for i in range(0,tamanho_t_filtro):
    plt.plot(x_MEG_SF_All_L[i],gamma_MX_H2O_MEG_L[i], '.')
for i in range(0, tamanho_t_filtro):
    plt.plot(x_MEG_SF_oti, gamma_MX_H2O_MEG_otimizacao[:,i], label=T_filtro[i])
plt.title("gamma_MX_H2O_MEG Calc")
plt.legend()

fig5 = plt.figure(5)
ax5 = fig5.add_subplot(111)
ax5.set_prop_cycle('color', colors)
plt.xlabel("x_MEG_SF")
plt.ylabel("Potencial")
for i in range(0,len(T_pred)):
    plt.plot(x_MEG, pot_calc_pred[:,i], label=T_pred[i])
plt.title("Potenciais Quimicos Pred")
plt.legend()

fig6 = plt.figure(6)
ax6 = fig6.add_subplot(111)
ax6.set_prop_cycle('color', colors)
plt.xlabel("x_MEG_SF")
plt.ylabel("b")
for i in range(0, len(T_pred)):
    plt.plot(x_MEG_SF_oti, b_MX_H2O_MEG_pred[:,i], label=T_pred[i])
plt.title("b_MX_H2O_MEG Pred")
plt.legend()

fig7 = plt.figure(7)
ax7 = fig7.add_subplot(111)
ax7.set_prop_cycle('color', colors)
plt.xlabel("x_MEG_SF")
plt.ylabel("gamma")
for i in range(0, len(T_pred)):
    plt.plot(x_MEG_SF_oti, gamma_MX_H2O_MEG_pred[:,i], label=T_pred[i])
plt.title("gamma_MX_H2O_MEG Pred")
plt.legend()

#plt.show()

from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages(f'saida_{nome_sal}_MI.pdf')
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

#Resultados
pd.ExcelWriter(f'saida_{nome_sal}_MI.xlsx', engine='openpyxl') 
#result = pd.read_excel('saida_sal_MultI.xlsx')

#Solubilidade puros
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

#Potencial quimico
x_MEG_df = pd.DataFrame(x_MEG_SF_All_tratado, columns=['x_MEG_SF'])
pot_df = pd.DataFrame(potencial_quimico_exc_tratado, columns=['Potencial exp'])
pot_calc_df = pd.DataFrame(pot_calc, columns=['Potencial calc'])
T_All_tratado_df = pd.DataFrame(T_All_tratado, columns=['T potencial'])
Thetas_pot_df = pd.DataFrame(Thetas_pot, columns=['Thetas pot'])
R_2_pot_df = pd.DataFrame({'R2 pot':[R_2_pot]})
AAD_pot_df = pd.DataFrame({'AAD pot':[AAD_pot]})
maxAD_pot_df = pd.DataFrame({'maxAD pot':[maxAD_pot]})
AARD_pot_df = pd.DataFrame({'AARD pot':[AARD_pot]})
maxARD_pot_df = pd.DataFrame({'maxARD pot':[maxARD_pot]})

#Solubilidade mistura + Gamma
x_MEG_SF_All_df = pd.DataFrame(x_MEG_SF_All, columns=['x_MEG_sol'])
w_MX_All_df = pd.DataFrame(w_MX_All, columns=['w MX'])
T_all_df2 = pd.DataFrame(T_All, columns=['T sol mix'])
b_MX_H2O_MEG_df = pd.DataFrame(b_MX_H2O_MEG, columns=['b_mix exp'])
gamma_MX_H2O_MEG_df = pd.DataFrame(gamma_MX_H2O_MEG, columns=['gamma exp'])
x_MEG_SF_oti_df = pd.DataFrame(x_MEG_SF_oti, columns=['x MEG SF calc'])
w_MEG_SF_oti_df = pd.DataFrame(w_MEG_SF_oti, columns=['w MEG SF calc'])
b_MX_H2O_MEG_otimizacao_df = pd.DataFrame(b_MX_H2O_MEG_otimizacao)
gamma_MX_H2O_MEG_otimizacao_df = pd.DataFrame(gamma_MX_H2O_MEG_otimizacao)
T_filtro_df = pd.DataFrame(T_filtro, columns=['Temperaturas isotermas'])
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


result = pd.concat([b_bin_H2O_df,T_bin_H2O_df, b_bin_MEG_df, T_bin_MEG_df,Thetas_H2O_df,Thetas_MEG_df,r2_H2O_df,r2_MEG_df,x_MEG_df,pot_df, pot_calc_df, T_All_tratado_df, Thetas_pot_df, R_2_pot_df, AAD_pot_df, maxAD_pot_df, AARD_pot_df, maxARD_pot_df, Beta_1_MX_MEG_df, Beta_1_MX_H2O_MEG_df, T_all_df1, x_MEG_SF_All_df, w_MX_All_df,T_all_df2, b_MX_H2O_MEG_df,b_MX_estimado_df, gamma_MX_H2O_MEG_df, gamma_MX_estimado_df, R_2_b_df, AAD_b_df, maxAD_b_df, AARD_b_df, maxARD_b_df, R_2_g_df, AAD_g_df,  maxAD_g_df, AARD_g_df, maxARD_g_df, x_MEG_SF_oti_df, w_MEG_SF_oti_df, T_filtro_df, b_MX_H2O_MEG_otimizacao_df, gamma_MX_H2O_MEG_otimizacao_df], axis=1)
result.to_excel(f'saida_{nome_sal}_MI.xlsx', sheet_name = nome_sal, header=True, index=False)


# %%
