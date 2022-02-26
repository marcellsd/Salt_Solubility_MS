#%%
import sys
from tkinter.tix import COLUMN

sys.path.append('Modulos')

import pandas as pd
from numpy import array, linspace, log, isnan, empty, exp, mean, concatenate, atleast_2d
import matplotlib.pyplot as plt
from tkinter import filedialog, Tk

from Modulos import transformacoes
from Modulos import pitzer
from Modulos import solubilidade
from Modulos import excesso
from Modulos import estatistica
from Modulos import filtro_dados
from matplotlib.backends.backend_pdf import PdfPages



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
w_MX_H2O_bin_L, w_H2O_bin_L = [], []

w_MEG_bin_L, w_MX_MEG_bin_L = [], []

w_H2O_tern_L, w_MEG_tern_L, w_MX_tern_L, w_H2O_SF_tern_L, w_MEG_SF_tern_L = [], [], [], [], []

w_H2O_All_L, w_MEG_All_L, w_MX_All_L, w_H2O_SF_All_L, w_MEG_SF_All_L, T_All_L = [], [], [], [], [], []   

w_H2O_SF_in_bin_H2O_L, w_H2O_SF_in_bin_MEG_L = [], []
w_MEG_SF_in_bin_MEG_L, w_MEG_SF_in_bin_H2O_L = [], []


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
        w_H2O_SF_in_bin_H2O_L.append(w_H2O_SF[i])
        w_MEG_SF_in_bin_H2O_L.append(1-w_H2O_SF[i])
    
    elif w_H2O_SF[i] <= 0.01:
        w_MEG_bin_L.append((1-w_H2O_SF[i])*(1 - w_MX_Todos[i]))
        w_MX_MEG_bin_L.append(w_MX_Todos[i])
        w_MEG_SF_in_bin_MEG_L.append(1-w_H2O_SF[i])
        w_H2O_SF_in_bin_MEG_L.append(w_H2O_SF[i])
    else:
        w_H2O_tern_L.append(w_H2O_SF[i]*(1 - w_MX_Todos[i]))
        w_MX_tern_L.append(w_MX_Todos[i])
        w_MEG_tern_L.append(1 - w_H2O_SF[i]*(1-w_MX_Todos[i]) - w_MX_Todos[i])
        w_H2O_SF_tern_L.append(w_H2O_SF[i])
        w_MEG_SF_tern_L.append(1-w_H2O_SF[i])  

#Substituir valores 0

for i in range (0, len(w_MX_MEG_bin_L)):
    if w_MX_MEG_bin_L[i] == 0:
        w_MX_MEG_bin_L[i] = 1e-6

#Transformando de lista para array
w_MX_H2O_bin, w_H2O_bin  = array(w_MX_H2O_bin_L), array(w_H2O_bin_L)

w_MEG_bin, w_MX_MEG_bin = array(w_MEG_bin_L, dtype="object"), array(w_MX_MEG_bin_L, dtype="object")

w_H2O_tern, w_MEG_tern, w_MX_tern, w_H2O_SF_tern = array(w_H2O_tern_L), array(w_MEG_tern_L), array(w_MX_tern_L), array(w_H2O_SF_tern_L)

w_H2O_All, w_MEG_All, w_MX_All, w_H2O_SF_All, w_MEG_SF_All, T_All = array(w_H2O_All_L), array(w_MEG_All_L), array(w_MX_All_L), array(w_H2O_SF_All_L), array(w_MEG_SF_All_L), array(T_All_L)

w_MEG_SF_tern = array(w_MEG_SF_tern_L)

w_H2O_SF_in_bin_H2O, w_MEG_SF_in_bin_H2O, w_MEG_SF_in_bin_MEG, w_H2O_SF_in_bin_MEG = array(w_H2O_SF_in_bin_H2O_L), array(w_MEG_SF_in_bin_H2O_L), array(w_MEG_SF_in_bin_MEG_L), array(w_H2O_SF_in_bin_MEG_L)

#Calculo das médias
w_MX_H2O_bin_mean = array(mean(w_MX_H2O_bin))
w_MX_MEG_bin_mean = array(mean(w_MX_MEG_bin))
w_H2O_SF_in_bin_H2O_mean = array(mean(w_H2O_SF_in_bin_H2O))
w_MEG_SF_in_bin_H2O_mean = array(mean(w_MEG_SF_in_bin_H2O))
w_MEG_SF_in_bin_MEG_mean = array(mean(w_MEG_SF_in_bin_MEG))
w_H2O_SF_in_bin_MEG_mean = array(mean(w_H2O_SF_in_bin_MEG))

#Juntando arrays
w_H2O_SF_All_True = concatenate([atleast_2d(a) for a in [w_H2O_SF_in_bin_H2O_mean,w_H2O_SF_tern,w_H2O_SF_in_bin_MEG_mean]])


w_MEG_SF_All_True = concatenate([atleast_2d(a) for a in [w_MEG_SF_in_bin_H2O_mean,w_MEG_SF_tern,w_MEG_SF_in_bin_MEG_mean]])

w_MX_All_True = concatenate([atleast_2d(a) for a in [w_MX_H2O_bin_mean,w_MX_tern,w_MX_MEG_bin_mean]])

for i in range (0, len(w_H2O_SF_All_True)):
    if w_H2O_SF_All_True[i] >= 1e-3:
        w_H2O_SF_All_True[i] == 0
for i in range (0, len(w_MEG_SF_All_True)):
    if w_MEG_SF_All_True[i] >= 1e-3:
        w_MEG_SF_All_True[i] == 0

#Transformando Fração Mássica para Molar SF

x_MEG_SF_All = transformacoes.x_MEG_SF(w_H2O_SF_All_True)
x_H2O_SF_All = (1-x_MEG_SF_All)

x_MEG_SF_tern = transformacoes.x_MEG_SF(w_H2O_SF_tern)
x_H2O_SF_tern = (1-x_MEG_SF_tern)

#Calculo das molalidades

b_MX_H2O_bin = transformacoes.fracao_massica_para_molalidade(MM_MX, w_MX_H2O_bin)
b_MX_MEG_bin = transformacoes.fracao_massica_para_molalidade(MM_MX, w_MX_MEG_bin)
b_MX_H2O_mean = array(mean(b_MX_H2O_bin))
b_MX_MEG_mean = array(mean(b_MX_MEG_bin))


#Pitzer
A_phi_H2O = pitzer.A_phi_H2O(T_filtro)
rho_H2O, rho_MEG, epsilon_r_H2O, epsilon_r_MEG = pitzer.rho_and_epsilon_r_of_pures(T_filtro)
A_phi_MEG = pitzer.A_phi_MEG(A_phi_H2O, rho_H2O, epsilon_r_H2O, rho_MEG, epsilon_r_MEG)
rho_H2O_MEG, epsilon_r_H2O_MEG = pitzer.rho_and_epsilon_r_of_mixing(rho_H2O, rho_MEG, epsilon_r_H2O, epsilon_r_MEG, 
                                                                    w_H2O_SF_All_True, x_MEG_SF_All, x_H2O_SF_All, w_MEG_SF_All_True, T_filtro)

A_phi_H2O_MEG = pitzer.A_phi_of_Mixing(rho_H2O_MEG, rho_H2O, epsilon_r_H2O_MEG, epsilon_r_H2O, T_filtro)

Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O = pitzer.pitzerParameters_MX_H2O(Beta_0, Beta_1, Cphi, T_filtro, Zm, Rm, Rx)

ln_gamma_MX_H2O = pitzer.ln_gamma_MX_SI(nu_M, nu_X, Zm, Zx, b_MX_H2O_mean, A_phi_H2O, Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O)


Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG = pitzer.pitzerParameter_MX_MEG_SI(nu_M, nu_X, A_phi_MEG, Beta_0_MX_H2O, Beta_1_MX_H2O,
                                                                          Cphi_MX_H2O, ln_gamma_MX_H2O, b_MX_H2O_mean, b_MX_MEG_mean, 
                                                                          T_filtro, Zm, Zx, DeltaG_transf)



Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, Cphi_MX_H2O_MEG = pitzer.Pitzer_parameters_in_mixtures_SI(Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O, 
                                                                                             Beta_1_MX_MEG, epsilon_r_H2O, epsilon_r_MEG, epsilon_r_H2O_MEG)
#%%

#Cálculo do Potencial Químico de excesso


b_MX_H2O_MEG = transformacoes.fracao_massica_para_molalidade(MM_MX, w_MX_All_True)

ln_gamma_MX_MEG = pitzer.ln_gamma_MX_SI(nu_M, nu_X, Zm, Zx, b_MX_MEG_mean, A_phi_MEG, Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG)

ln_gamma_MX_H2O_MEG = pitzer.ln_gamma_MX_SI(nu_M, nu_X, Zm, Zx, b_MX_H2O_MEG, A_phi_H2O_MEG, Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, Cphi_MX_H2O_MEG, False)  
ln_gamma_MX_H2O_MEG_ideal = excesso.ln_gamma_MX_H2O_MEG_ideal(x_MEG_SF_All, ln_gamma_MX_H2O, ln_gamma_MX_MEG)
ln_gamma_MX_H2O_MEG_exc = excesso.ln_gamma_MX_H2O_MEG_exc(ln_gamma_MX_H2O_MEG_ideal, ln_gamma_MX_H2O_MEG)

ln_b_MX_H2O_MEG = log(b_MX_H2O_MEG)
ln_b_MX_H2O_MEG_ideal = excesso.ln_b_MX_H2O_MEG_ideal(x_MEG_SF_All ,b_MX_H2O_mean, b_MX_MEG_mean)
ln_b_MX_H2O_MEG_exc = excesso.ln_b_MX_H2O_MEG_exc(ln_b_MX_H2O_MEG_ideal, ln_b_MX_H2O_MEG)

potencial_quimico_exc = excesso.potencial_quimico_exc(ln_b_MX_H2O_MEG_exc, ln_gamma_MX_H2O_MEG_exc)

#Quantidade de parametros potencial quimico (escolher 3 ou 2)
n = 3
     
Thetas_pot = excesso.regressao_n_linear_potencial_quimico_SI(potencial_quimico_exc, x_MEG_SF_All, n)

pot_calc = excesso.func_potencial_quimico_exc_SI(x_MEG_SF_All, Thetas_pot, n)

x_MEG = linspace( 0.0,1.0, 100, endpoint=True)

pot_calc_plot = excesso.func_potencial_quimico_exc_SI(x_MEG, Thetas_pot, n)

#Estatistica potencial
R_2_pot = estatistica.coef_determinacao2(potencial_quimico_exc, pot_calc)
AAD_pot, maxAD_pot, AARD_pot, maxARD_pot = estatistica.DAM_DMR(potencial_quimico_exc, pot_calc)


#Calculo solubilidade Mistura
numero_de_pontos = 31

x_MEG_SF_oti = linspace(0.0,1.0,numero_de_pontos,endpoint=True)
x_H2O_SF_oti = 1-x_MEG_SF_oti
w_MEG_SF_oti = transformacoes.w_MEG_SF(x_MEG_SF_oti)
w_H2O_SF_oti = 1-w_MEG_SF_oti

b_MX_H2O_MEG_otimizacao = empty(numero_de_pontos)
gamma_MX_H2O_MEG_otimizacao = empty(numero_de_pontos)

b_MX_H2O_MEG_otimizacao, gamma_MX_H2O_MEG_otimizacao = excesso.otimizacao_b_MX_H2O_MEG_SI(x_MEG_SF_oti, w_H2O_SF_oti, x_H2O_SF_oti, w_MEG_SF_oti, T_filtro[0], 
                                                                                                     Thetas_pot, b_MX_H2O_mean, b_MX_MEG_mean, 
                                                                                                     Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O,
                                                                                                     Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG,
                                                                                                     Zm, nu_M, nu_X, Zx, n)

b_MX_H2O_MEG_est = empty(numero_de_pontos)
gamma_MX_H2O_MEG_est = empty(numero_de_pontos)
b_MX_H2O_MEG_est, gamma_MX_H2O_MEG_est = excesso.otimizacao_b_MX_H2O_MEG_SI(x_MEG_SF_All[:,0], w_H2O_SF_All_True[:,0], x_H2O_SF_All[:,0], w_MEG_SF_All_True[:,0], T_filtro[0], 
                                                                                                     Thetas_pot, b_MX_H2O_mean, b_MX_MEG_mean, 
                                                                                                     Beta_0_MX_H2O[0], Beta_1_MX_H2O[0], Cphi_MX_H2O[0],
                                                                                                     Beta_0_MX_MEG[0], Beta_1_MX_MEG, Cphi_MX_MEG[0],
                                                                                                     Zm, nu_M, nu_X, Zx, n)

gamma_MX_H2O_MEG = exp(ln_gamma_MX_H2O_MEG)

#Estatistica solubilidade e gamma
R_2_b=estatistica.coef_determinacao2(b_MX_H2O_MEG, b_MX_H2O_MEG_est)
R_2_gamma=estatistica.coef_determinacao2(gamma_MX_H2O_MEG, gamma_MX_H2O_MEG_est)
AAD_b, maxAD_b, AARD_b, maxARD_b = estatistica.DAM_DMR(b_MX_H2O_MEG, b_MX_H2O_MEG_est)
AAD_g, maxAD_g, AARD_g, maxARD_g = estatistica.DAM_DMR(gamma_MX_H2O_MEG, gamma_MX_H2O_MEG_est)

print (f"R_2_b: {R_2_b:.4f}, R_2_gamma: {R_2_gamma:.4f}, R_Pot: {R_2_pot:.4f}")

fig1 = plt.figure(1)      
plt.figure(1)
plt.xlabel("x_MEG_SF")
plt.ylabel("Potencial")
plt.plot(x_MEG_SF_All, potencial_quimico_exc, '.')
plt.plot(x_MEG, pot_calc_plot, label=T_filtro)
plt.title("Potenciais Quimicos Calc")
plt.legend()


#Grafico b
fig2 = plt.figure(2)
plt.xlabel("x_MEG_SF")
plt.ylabel("b")
plt.plot(x_MEG_SF_All, b_MX_H2O_MEG, '.')
plt.plot(x_MEG_SF_oti, b_MX_H2O_MEG_otimizacao, label=T_filtro)
plt.title("b_MX_H2O_MEG Calc")
plt.legend()

#Grafico gamma
fig3 = plt.figure(3)
plt.xlabel("x_MEG_SF")
plt.ylabel("gamma")
plt.plot(x_MEG_SF_All,gamma_MX_H2O_MEG, '.')
plt.plot(x_MEG_SF_oti, gamma_MX_H2O_MEG_otimizacao, label=T_filtro)
plt.title("gamma_MX_H2O_MEG Calc")
plt.legend()

fig4 = plt.figure(4)
plt.xlabel("x_MEG_SF")
plt.ylabel("ln_gamma")
plt.plot(x_MEG_SF_All,ln_gamma_MX_H2O_MEG, '.')
plt.plot(x_MEG_SF_oti, log(gamma_MX_H2O_MEG_otimizacao), label=T_filtro)
plt.title("gamma_MX_H2O_MEG Calc")
plt.legend()

fig5 = plt.figure(5)
plt.xlabel("x_MEG_SF")
plt.ylabel("ln_b")
plt.plot(x_MEG_SF_All, ln_b_MX_H2O_MEG, '.')
plt.plot(x_MEG_SF_oti, log(b_MX_H2O_MEG_otimizacao), label=T_filtro)
plt.title("b_MX_H2O_MEG Calc")
plt.legend()

pp = PdfPages(f'saida_{nome_sal}_{int(T_filtro)}_SI.pdf')
plt.savefig(pp, format='pdf')
pp.savefig(fig1)
pp.savefig(fig2)
pp.savefig(fig3)
pp.savefig(fig4)
pp.savefig(fig5)
pp.close()

#Resultados
pd.ExcelWriter(f'saida_{nome_sal}_{int(T_filtro)}_SI.xlsx', engine='openpyxl') 
#result = pd.read_excel(f'{nome_sal}_results.xlsx')

#Pitzer

T_filtro_df = pd.DataFrame(T_filtro, columns=['T potencial'])
Beta_1_MX_MEG_df = pd.DataFrame(Beta_1_MX_MEG, columns=['Beta 1 MX MEG'])
Beta_1_MX_H2O_MEG_df = pd.DataFrame(Beta_1_MX_H2O_MEG, columns=['Beta 1 MX H2O MEG'])

#Potencial quimico
x_MEG_df = pd.DataFrame(x_MEG_SF_All, columns=['x_MEG_SF'])
pot_df = pd.DataFrame(potencial_quimico_exc, columns=['Potencial exp'])
pot_calc_df = pd.DataFrame(pot_calc, columns=['Potencial calc'])
Thetas_pot_df = pd.DataFrame(Thetas_pot, columns=['Thetas pot'])
R_2_pot_df = pd.DataFrame({'R2 pot':[R_2_pot]})
AAD_pot_df = pd.DataFrame({'AAD pot':[AAD_pot]})
maxAD_pot_df = pd.DataFrame({'maxAD pot':[maxAD_pot]})
AARD_pot_df = pd.DataFrame({'AARD pot':[AARD_pot]})
maxARD_pot_df = pd.DataFrame({'maxARD pot':[maxARD_pot]})

#Solubilidade mistura + Gamma
x_MEG_SF_All_df = pd.DataFrame(x_MEG_SF_All, columns=['x_MEG_sol'])
w_MX_All_df = pd.DataFrame(w_MX_All_True, columns=['w MX'])
b_MX_H2O_MEG_df = pd.DataFrame(b_MX_H2O_MEG, columns=['b_mix exp'])
gamma_MX_H2O_MEG_df = pd.DataFrame(gamma_MX_H2O_MEG, columns=['gamma exp'])
x_MEG_SF_oti_df = pd.DataFrame(x_MEG_SF_oti, columns=['x MEG SF calc'])
w_MEG_SF_oti_df = pd.DataFrame(w_MEG_SF_oti, columns=['w MEG SF calc'])
b_MX_H2O_MEG_otimizacao_df = pd.DataFrame(b_MX_H2O_MEG_otimizacao, columns=['b_MX_H2O_MEG_otimizado'])
gamma_MX_H2O_MEG_otimizacao_df = pd.DataFrame(gamma_MX_H2O_MEG_otimizacao, columns=['gamma_MX_H2O_MEG_otimizado'])
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

result = pd.concat([T_filtro_df,x_MEG_df,pot_df, pot_calc_df, Thetas_pot_df, R_2_pot_df, AAD_pot_df, maxAD_pot_df, AARD_pot_df, maxARD_pot_df, Beta_1_MX_MEG_df, Beta_1_MX_H2O_MEG_df, x_MEG_SF_All_df, w_MX_All_df, b_MX_H2O_MEG_df,b_MX_estimado_df, gamma_MX_H2O_MEG_df, gamma_MX_estimado_df, R_2_b_df, AAD_b_df, maxAD_b_df, AARD_b_df, maxARD_b_df, R_2_g_df, AAD_g_df,  maxAD_g_df, AARD_g_df, maxARD_g_df, x_MEG_SF_oti_df, w_MEG_SF_oti_df, b_MX_H2O_MEG_otimizacao_df, gamma_MX_H2O_MEG_otimizacao_df], axis=1)
result.to_excel(f'saida_{nome_sal}_{int(T_filtro)}_SI.xlsx', sheet_name = nome_sal, header=True, index=False)

#plt.show()
# %%

