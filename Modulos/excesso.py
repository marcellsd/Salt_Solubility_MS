# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 10:47:40 2021

@author: marcell
"""
from numpy import log, empty, var, exp
from pyswarm import pso
from scipy.optimize import minimize, fsolve
import solubilidade
import pitzer

def ln_b_MX_H2O_MEG_ideal(x_MEG_SF_tern, b_MX_H2O, b_MX_MEG):
    try:
        tamanho_v = len(x_MEG_SF_tern)
        tamanho_b = len(b_MX_H2O)
        ln_b_MX_H2O_MEG_ideal = empty(tamanho_v)
        if tamanho_v == 1 and tamanho_b == 1:
            ln_b_MX_H2O_MEG_ideal = (1-x_MEG_SF_tern)*log(b_MX_H2O)+x_MEG_SF_tern*log(b_MX_MEG)
        elif tamanho_b == 1:
            for i in range(0,tamanho_v):
                ln_b_MX_H2O_MEG_ideal[i] = (1-x_MEG_SF_tern[i])*log(b_MX_H2O[0])+x_MEG_SF_tern[i]*log(b_MX_MEG[0])
        else:
            for i in range(0,tamanho_v):
                ln_b_MX_H2O_MEG_ideal[i] = (1-x_MEG_SF_tern[i])*log(b_MX_H2O[i])+x_MEG_SF_tern[i]*log(b_MX_MEG[i])
    except:
        tamanho_v = len(x_MEG_SF_tern)
        ln_b_MX_H2O_MEG_ideal = empty(tamanho_v)
        for i in range(0,tamanho_v):
                ln_b_MX_H2O_MEG_ideal[i] = (1-x_MEG_SF_tern[i])*log(b_MX_H2O)+x_MEG_SF_tern[i]*log(b_MX_MEG)
    return ln_b_MX_H2O_MEG_ideal

def ln_gamma_MX_H2O_MEG_ideal(x_MEG_SF_tern, ln_gamma_MX_H2O, ln_gamma_MX_MEG):
    try:
        tamanho_v = len(x_MEG_SF_tern)
        tamanho_g = len(ln_gamma_MX_H2O)
        ln_gamma_MX_H2O_MEG_ideal = empty(tamanho_v)
        if tamanho_v == 1 and tamanho_g == 1:
            ln_gamma_MX_H2O_MEG_ideal = (1-x_MEG_SF_tern)*ln_gamma_MX_H2O + x_MEG_SF_tern*ln_gamma_MX_MEG
        elif tamanho_g == 1:
            for i in range(0,tamanho_v):
                ln_gamma_MX_H2O_MEG_ideal[i] = (1-x_MEG_SF_tern[i])*ln_gamma_MX_H2O[0] + x_MEG_SF_tern[i] * ln_gamma_MX_MEG[0]
        else:
            for i in range(0,tamanho_v):
                ln_gamma_MX_H2O_MEG_ideal[i] = (1-x_MEG_SF_tern[i])*ln_gamma_MX_H2O[i] + x_MEG_SF_tern[i]*ln_gamma_MX_MEG[i]
    except:
        tamanho_v = len(x_MEG_SF_tern)
        ln_gamma_MX_H2O_MEG_ideal = empty(tamanho_v)
        for i in range(0,tamanho_v):
                ln_gamma_MX_H2O_MEG_ideal[i] = (1-x_MEG_SF_tern[i])*ln_gamma_MX_H2O + x_MEG_SF_tern[i] * ln_gamma_MX_MEG
    return ln_gamma_MX_H2O_MEG_ideal

def ln_b_MX_H2O_MEG_exc(ln_b_MX_H2O_MEG_ideal, ln_b_MX_H2O_MEG):
    tamanho_v = len(ln_b_MX_H2O_MEG_ideal)
    ln_b_MX_H2O_MEG_exc = empty(tamanho_v)
    if tamanho_v == 1:
        ln_b_MX_H2O_MEG_exc = ln_b_MX_H2O_MEG - ln_b_MX_H2O_MEG_ideal
    else:
        for i in range(0,tamanho_v):
            ln_b_MX_H2O_MEG_exc[i] = ln_b_MX_H2O_MEG[i] - ln_b_MX_H2O_MEG_ideal[i]
    return ln_b_MX_H2O_MEG_exc

def ln_gamma_MX_H2O_MEG_exc(ln_gamma_MX_H2O_MEG_ideal, ln_gamma_MX_H2O_MEG):
    tamanho_v = len(ln_gamma_MX_H2O_MEG_ideal)
    ln_gamma_MX_H2O_MEG_exc = empty(tamanho_v)
    if tamanho_v == 1:
         ln_gamma_MX_H2O_MEG_exc =  ln_gamma_MX_H2O_MEG - ln_gamma_MX_H2O_MEG_ideal
    else:
        for i in range (0, tamanho_v):
           ln_gamma_MX_H2O_MEG_exc[i] =  ln_gamma_MX_H2O_MEG[i] - ln_gamma_MX_H2O_MEG_ideal[i]
    return ln_gamma_MX_H2O_MEG_exc

def potencial_quimico_exc(ln_b_MX_H2O_MEG_exc, ln_gamma_MX_H2O_MEG_exc):
    tamanho_v = len(ln_b_MX_H2O_MEG_exc)
    potencial_quimico_exc = empty(tamanho_v)
    
    if tamanho_v == 1:
        potencial_quimico_exc = ln_b_MX_H2O_MEG_exc + ln_gamma_MX_H2O_MEG_exc
    else:
        for i in range(0, tamanho_v):
            potencial_quimico_exc[i] = ln_b_MX_H2O_MEG_exc[i] + ln_gamma_MX_H2O_MEG_exc[i]
    return potencial_quimico_exc

def func_potencial_quimico_exc(T, x_MEG_SF_tern, Thetas_pot, n):
    if n == 4:
        pot_exc = x_MEG_SF_tern*(1-x_MEG_SF_tern)*(Thetas_pot[0]+Thetas_pot[1]*T+(Thetas_pot[2]+Thetas_pot[3]*T)*x_MEG_SF_tern)
    elif n == 6:                                               
        pot_exc = x_MEG_SF_tern*(1-x_MEG_SF_tern)*(Thetas_pot[0]+Thetas_pot[1]*T+(Thetas_pot[2]+Thetas_pot[3]*T)*x_MEG_SF_tern+(Thetas_pot[4]+Thetas_pot[5]*T)*x_MEG_SF_tern**2)
    return pot_exc

def fobjetivo(Thetas_pot, potencial_quimico_exc, T_tern, x_MEG_SF_tern, n):
    f = func_potencial_quimico_exc(T_tern, x_MEG_SF_tern, Thetas_pot, n)
    sol = 0.0
    TF = len(f)
    for i in range(0,TF):
        sol = sol + (potencial_quimico_exc[i] - f[i])**2
    return sol     
  

def regressao_n_linear_potencial_quimico(potencial_quimico_exc, T_tern, x_MEG_SF_tern, n):
    if n ==6:
        lb = [-100.0, -100.0, -100.0, -100.0, -100.0, -100.0]
        ub = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
    elif n == 4:
        lb = [-100.0, -100.0, -100.0, -100.0]
        ub = [100.0, 100.0, 100.0, 100.0]
    
    #Otimizador Heurístico
    param, fopt = pso(fobjetivo, lb, ub, ieqcons=[], f_ieqcons=None, args=(potencial_quimico_exc, T_tern, x_MEG_SF_tern, n), kwargs={}, swarmsize=100,
                 omega=0.5, phip=0.8, phig=0.8, maxiter=5000, minstep=1e-6, minfunc=1e-6, debug=False)
    
    x0 = (param)

    #Otimizador Determinístico
    Thetas_pot = minimize(fobjetivo,x0,args=(potencial_quimico_exc, T_tern, x_MEG_SF_tern, n),method='Nelder-Mead', tol = 1e-7, options = {'xatol':1e-7, 'disp':True, 'maxiter': 2000})
    return Thetas_pot.x
    
def otimizacao_b_MX_H2O_MEG(x_MEG_SF, w_H2O_SF, x_H2O_SF, w_MEG_SF, T, 
                            Thetas_pot, Thetas_H2O, Thetas_MEG, 
                            Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O,
                            Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG,
                            Zm, nu_M, nu_X, Zx, n):
    
    pot_exc_calc = func_potencial_quimico_exc(T, x_MEG_SF, Thetas_pot, n)
    b_MX_H2O = solubilidade.func_solubilidade(T,Thetas_H2O)
    b_MX_MEG = solubilidade.func_solubilidade(T,Thetas_MEG)
    
    ln_b_MX_H2O_MEG_ideal_calc = ln_b_MX_H2O_MEG_ideal(x_MEG_SF, b_MX_H2O, b_MX_MEG)
    
    #Pre Pitzer
    A_phi_H2O = pitzer.A_phi_H2O(T)
    rho_H2O, rho_MEG, epsilon_r_H2O, epsilon_r_MEG = pitzer.rho_and_epsilon_r_of_pures(T)
    A_phi_MEG = pitzer.A_phi_MEG(A_phi_H2O, rho_H2O, epsilon_r_H2O, rho_MEG, epsilon_r_MEG)
    ln_gamma_MX_H2O = pitzer.ln_gamma_MX(nu_M, nu_X, Zm, Zx, b_MX_H2O, A_phi_H2O, Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O)
    ln_gamma_MX_MEG = pitzer.ln_gamma_MX(nu_M, nu_X, Zm, Zx, b_MX_MEG, A_phi_MEG, Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG)
    
    #Pitzer Ideal
    ln_gamma_MX_H2O_MEG_ideal_calc = ln_gamma_MX_H2O_MEG_ideal(x_MEG_SF, ln_gamma_MX_H2O, ln_gamma_MX_MEG)
    
    #Pitzer Calculado
    rho_H2O_MEG, epsilon_r_H2O_MEG = pitzer.rho_and_epsilon_r_of_mixing(rho_H2O, rho_MEG, epsilon_r_H2O, epsilon_r_MEG, w_H2O_SF, x_MEG_SF, x_H2O_SF, w_MEG_SF, T)
    A_phi_H2O_MEG = pitzer.A_phi_of_Mixing(rho_H2O_MEG, rho_H2O, epsilon_r_H2O_MEG, epsilon_r_H2O, T)
    Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, Cphi_MX_H2O_MEG = pitzer.Pitzer_parameters_in_mixtures_exc(Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O, Beta_1_MX_MEG, epsilon_r_H2O, epsilon_r_MEG, epsilon_r_H2O_MEG, x_MEG_SF)
    
    def minimizacao(b_MX_H2O_MEG):
        ln_gamma_MX_H2O_MEG = pitzer.ln_gamma_MX(nu_M, nu_X, Zm, Zx, b_MX_H2O_MEG, A_phi_H2O_MEG, Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, Cphi_MX_H2O_MEG)
        ln_b_MX_H2O_MEG = log(b_MX_H2O_MEG)
        return (pot_exc_calc - (ln_b_MX_H2O_MEG - ln_b_MX_H2O_MEG_ideal_calc) - (ln_gamma_MX_H2O_MEG - ln_gamma_MX_H2O_MEG_ideal_calc))**2
        
    x0 = empty(len(ln_b_MX_H2O_MEG_ideal_calc))
    for i in range(0, len(ln_b_MX_H2O_MEG_ideal_calc)):
        x0[i] = exp(ln_b_MX_H2O_MEG_ideal_calc[i])
    
    b_MX_H2O_MEG_otimizado = fsolve(minimizacao,x0)
    
    ln_gamma_MX_H2O_MEG_otimizado = pitzer.ln_gamma_MX(nu_M, nu_X, Zm, Zx, b_MX_H2O_MEG_otimizado, A_phi_H2O_MEG, Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, Cphi_MX_H2O_MEG)
    gamma_MX_H2O_MEG_otimizado = exp(ln_gamma_MX_H2O_MEG_otimizado)
    return b_MX_H2O_MEG_otimizado, gamma_MX_H2O_MEG_otimizado
    

def otimizacao_b_MX_H2O_MEG_estatistica(x_MEG_SF, w_H2O_SF, x_H2O_SF, w_MEG_SF, T, 
                        Thetas_pot, Thetas_H2O, Thetas_MEG, 
                        Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O,
                        Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG,
                        Zm, nu_M, nu_X, Zx, n):

    pot_exc_calc = func_potencial_quimico_exc(T, x_MEG_SF, Thetas_pot, n)
    b_MX_H2O = solubilidade.func_solubilidade(T,Thetas_H2O)
    b_MX_MEG = solubilidade.func_solubilidade(T,Thetas_MEG)
    
    ln_b_MX_H2O_MEG_ideal_calc = ln_b_MX_H2O_MEG_ideal(x_MEG_SF, b_MX_H2O, b_MX_MEG)
    
    #Pre Pitzer
    A_phi_H2O = pitzer.A_phi_H2O(T)
    rho_H2O, rho_MEG, epsilon_r_H2O, epsilon_r_MEG = pitzer.rho_and_epsilon_r_of_pures(T)
    A_phi_MEG = pitzer.A_phi_MEG(A_phi_H2O, rho_H2O, epsilon_r_H2O, rho_MEG, epsilon_r_MEG)
    ln_gamma_MX_H2O = pitzer.ln_gamma_MX(nu_M, nu_X, Zm, Zx, b_MX_H2O, A_phi_H2O, Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O)
    ln_gamma_MX_MEG = pitzer.ln_gamma_MX(nu_M, nu_X, Zm, Zx, b_MX_MEG, A_phi_MEG, Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG)
    
    #Pitzer Ideal
    ln_gamma_MX_H2O_MEG_ideal_calc = ln_gamma_MX_H2O_MEG_ideal(x_MEG_SF, ln_gamma_MX_H2O, ln_gamma_MX_MEG)
    
    #Pitzer Calculado
    rho_H2O_MEG, epsilon_r_H2O_MEG = pitzer.rho_and_epsilon_r_of_mixing(rho_H2O, rho_MEG, epsilon_r_H2O, epsilon_r_MEG, w_H2O_SF, x_MEG_SF, x_H2O_SF, w_MEG_SF, T)
    A_phi_H2O_MEG = pitzer.A_phi_of_Mixing(rho_H2O_MEG, rho_H2O, epsilon_r_H2O_MEG, epsilon_r_H2O, T)
    Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, Cphi_MX_H2O_MEG = pitzer.Pitzer_parameters_in_mixtures(Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O, Beta_1_MX_MEG, epsilon_r_H2O, epsilon_r_MEG, epsilon_r_H2O_MEG)
    
    def minimizacao(b_MX_H2O_MEG):
        ln_gamma_MX_H2O_MEG = pitzer.ln_gamma_MX(nu_M, nu_X, Zm, Zx, b_MX_H2O_MEG, A_phi_H2O_MEG, Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, Cphi_MX_H2O_MEG)
        ln_b_MX_H2O_MEG = log(b_MX_H2O_MEG)
        resp = empty(len(ln_gamma_MX_H2O_MEG))
        for i in range(len(resp)):
            resp[i] = (pot_exc_calc[i] - (ln_b_MX_H2O_MEG[i] - ln_b_MX_H2O_MEG_ideal_calc[i]) - (ln_gamma_MX_H2O_MEG[i] - ln_gamma_MX_H2O_MEG_ideal_calc[i]))**2
        return resp
        
    x0 = exp(ln_b_MX_H2O_MEG_ideal_calc)
    b_MX_H2O_MEG_otimizado = fsolve(minimizacao,x0)
    
    ln_gamma_MX_H2O_MEG_otimizado = pitzer.ln_gamma_MX(nu_M, nu_X, Zm, Zx, b_MX_H2O_MEG_otimizado, A_phi_H2O_MEG, Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, Cphi_MX_H2O_MEG)
    gamma_MX_H2O_MEG_otimizado = exp(ln_gamma_MX_H2O_MEG_otimizado)
    return b_MX_H2O_MEG_otimizado, gamma_MX_H2O_MEG_otimizado





def ln_gamma_MX_H2O_MEG_exc_SI(ln_gamma_MX_H2O_MEG_ideal, ln_gamma_MX_H2O_MEG):
    tamanho_v = len(ln_gamma_MX_H2O_MEG_ideal)
    ln_gamma_MX_H2O_MEG_exc = empty(tamanho_v)
    if tamanho_v == 1:
         ln_gamma_MX_H2O_MEG_exc =  ln_gamma_MX_H2O_MEG - ln_gamma_MX_H2O_MEG_ideal
    else:
        for i in range (0, tamanho_v):
           ln_gamma_MX_H2O_MEG_exc[i] =  ln_gamma_MX_H2O_MEG[i] - ln_gamma_MX_H2O_MEG_ideal[i]
    return ln_gamma_MX_H2O_MEG_exc



def func_potencial_quimico_exc_SI(x_MEG_SF, Thetas_pot, n):
    if n == 2:
        pot_exc = x_MEG_SF*(1-x_MEG_SF)*(Thetas_pot[0]+(Thetas_pot[1])*x_MEG_SF)
    elif n == 3:                                               
        pot_exc = x_MEG_SF*(1-x_MEG_SF)*(Thetas_pot[0]+(Thetas_pot[1])*x_MEG_SF+(Thetas_pot[2])*x_MEG_SF**2)
    return pot_exc

def fobjetivo_SI(Thetas_pot, potencial_quimico_exc, x_MEG_SF, n):
    f = func_potencial_quimico_exc_SI(x_MEG_SF, Thetas_pot, n)
    sol = 0.0
    TF = len(f)
    for i in range(0,TF):
        sol = sol + (potencial_quimico_exc[i] - f[i])**2
    return sol  


def potencial_quimico_exc_SI(ln_b_MX_H2O_MEG_exc, ln_gamma_MX_H2O_MEG_exc):
    tamanho_v = len(ln_b_MX_H2O_MEG_exc)
    potencial_quimico_exc = empty(tamanho_v)
    
    if tamanho_v == 1:
        potencial_quimico_exc = ln_b_MX_H2O_MEG_exc + ln_gamma_MX_H2O_MEG_exc
    else:
        for i in range(0, tamanho_v):
            potencial_quimico_exc[i] = ln_b_MX_H2O_MEG_exc[i] + ln_gamma_MX_H2O_MEG_exc[i]
    return potencial_quimico_exc


def regressao_n_linear_potencial_quimico_SI(potencial_quimico_exc, x_MEG_SF, n):
    if n ==3:
        lb = [-100.0, -100.0, -100.0]
        ub = [100.0, 100.0, 100.0]
    elif n == 2:
        lb = [-100.0, -100.00]
        ub = [100.0, 100.0]
    
    #Otimizador Heurístico
    param, fopt = pso(fobjetivo_SI, lb, ub, ieqcons=[], f_ieqcons=None, args=(potencial_quimico_exc, x_MEG_SF, n), kwargs={}, swarmsize=100,
                 omega=0.5, phip=0.8, phig=0.8, maxiter=5000, minstep=1e-6, minfunc=1e-6, debug=False)
    
    x0 = (param)

    #Otimizador Determinístico
    Thetas_pot = minimize(fobjetivo_SI,x0,args=(potencial_quimico_exc, x_MEG_SF, n),method='Nelder-Mead', tol = 1e-7, options = {'xatol':1e-7, 'disp':True, 'maxiter': 2000})
    return Thetas_pot.x

def otimizacao_b_MX_H2O_MEG_SI(x_MEG_SF, w_H2O_SF, x_H2O_SF, w_MEG_SF, T, 
                            Thetas_pot, b_MX_H2O, b_MX_MEG, 
                            Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O,
                            Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG,
                            Zm, nu_M, nu_X, Zx, n):
    
    pot_exc_calc = func_potencial_quimico_exc_SI(x_MEG_SF, Thetas_pot, n)
    
    ln_b_MX_H2O_MEG_ideal_calc = ln_b_MX_H2O_MEG_ideal(x_MEG_SF, b_MX_H2O, b_MX_MEG)
    
    #Pre Pitzer
    A_phi_H2O = pitzer.A_phi_H2O(T)
    rho_H2O, rho_MEG, epsilon_r_H2O, epsilon_r_MEG = pitzer.rho_and_epsilon_r_of_pures(T)
    A_phi_MEG = pitzer.A_phi_MEG(A_phi_H2O, rho_H2O, epsilon_r_H2O, rho_MEG, epsilon_r_MEG)
    ln_gamma_MX_H2O = pitzer.ln_gamma_MX_SI(nu_M, nu_X, Zm, Zx, b_MX_H2O, A_phi_H2O, Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O)
    ln_gamma_MX_MEG = pitzer.ln_gamma_MX_SI(nu_M, nu_X, Zm, Zx, b_MX_MEG, A_phi_MEG, Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG)
    
    #Pitzer Ideal
    ln_gamma_MX_H2O_MEG_ideal_calc = ln_gamma_MX_H2O_MEG_ideal(x_MEG_SF, ln_gamma_MX_H2O, ln_gamma_MX_MEG)
    
    #Pitzer Calculado
    rho_H2O_MEG, epsilon_r_H2O_MEG = pitzer.rho_and_epsilon_r_of_mixing(rho_H2O, rho_MEG, epsilon_r_H2O, epsilon_r_MEG, w_H2O_SF, x_MEG_SF, x_H2O_SF, w_MEG_SF, T)
    A_phi_H2O_MEG = pitzer.A_phi_of_Mixing(rho_H2O_MEG, rho_H2O, epsilon_r_H2O_MEG, epsilon_r_H2O, T)
    Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, Cphi_MX_H2O_MEG = pitzer.Pitzer_parameters_in_mixtures_exc_SI(Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O, Beta_1_MX_MEG, epsilon_r_H2O, epsilon_r_MEG, epsilon_r_H2O_MEG, x_MEG_SF)
    
    def minimizacao(b_MX_H2O_MEG):
        ln_gamma_MX_H2O_MEG = pitzer.ln_gamma_MX_SI(nu_M, nu_X, Zm, Zx, b_MX_H2O_MEG, A_phi_H2O_MEG, Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, Cphi_MX_H2O_MEG, False)
        ln_b_MX_H2O_MEG = log(b_MX_H2O_MEG)
        min = (pot_exc_calc - (ln_b_MX_H2O_MEG - ln_b_MX_H2O_MEG_ideal_calc) - (ln_gamma_MX_H2O_MEG - ln_gamma_MX_H2O_MEG_ideal_calc))**2
        return min


    x0 = empty(len(ln_b_MX_H2O_MEG_ideal_calc))
    for i in range(0, len(ln_b_MX_H2O_MEG_ideal_calc)):
        x0[i] = exp(ln_b_MX_H2O_MEG_ideal_calc[i])

    b_MX_H2O_MEG_otimizado = fsolve(minimizacao,x0)
    
    ln_gamma_MX_H2O_MEG_otimizado = pitzer.ln_gamma_MX_SI(nu_M, nu_X, Zm, Zx, b_MX_H2O_MEG_otimizado, A_phi_H2O_MEG, Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, Cphi_MX_H2O_MEG,False)
    gamma_MX_H2O_MEG_otimizado = exp(ln_gamma_MX_H2O_MEG_otimizado)
    return b_MX_H2O_MEG_otimizado, gamma_MX_H2O_MEG_otimizado