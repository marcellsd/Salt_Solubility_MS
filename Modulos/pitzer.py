# -*- coding: utf-8 -*-
"""
Created on Tue May 25 18:07:04 2021

@author: marcell
"""
from numpy import exp, log, sqrt, empty
from pyswarm import pso
from scipy.optimize import minimize, fsolve
from Modulos import constants

def rho_and_epsilon_r_of_pures(T):
    # Calcula as Densidades Absolutas e as Constantes Dielétricas
    # dos solventes puros
    # GREEN e PERRY, 2007
    
    tau = 1 - T/647.096
    rho_H2O = (17863 + 58606*tau**.35 - 95396*tau**(2/3) + 213890*tau - 141260*tau**(4/3)) * constants.MM_H2O
    rho_MEG = (1315 / (.25125**(1 + (1 - T/720)**.21868))) * constants.MM_MEG
    # ÅKERLÖF, 1932
    epsilon_r_H2O = 10**(1.9051 - .00205*(T-293.15))
    epsilon_r_MEG = 10**(1.5872 - .00224*(T-293.15))
    return rho_H2O, rho_MEG, epsilon_r_H2O, epsilon_r_MEG

def rho_and_epsilon_r_of_mixing(rho_H2O, rho_MEG, epsilon_r_H2O, epsilon_r_MEG, w_H2O_SF, x_MEG_SF, x_H2O_SF, w_MEG_SF, T):
    # Calcula a Densidade Absoluta e a Constante Dielétrica da
    # Regra de Mistura
    rho_H2O_MEG = 1/((w_H2O_SF/rho_H2O)+(w_MEG_SF/rho_MEG))
    # JOUYBAN, SOLTANPOUR e CHAN, 2004
    epsilon_r_H2O_MEG = exp(x_H2O_SF*log(epsilon_r_H2O) + x_MEG_SF*log(epsilon_r_MEG) +
        x_H2O_SF*x_MEG_SF/T * (153.6 + 57.3*(x_H2O_SF-x_MEG_SF)))
    return rho_H2O_MEG, epsilon_r_H2O_MEG

def A_phi_H2O(T):
    # Calcula a Constante de Debye-Hückel para a água
    # => CHEN, BRITT, et al., 1982
    return (-61.44534*exp((T-273.15)/273.15) + 2.864468*exp(2*(T-273.15)/273.15) + 183.5379*
            log(T/273.15) - .6820223*(T-273.15) + .0007875695*(T**2 - 273.15**2) +
            58.95788*(273.15/T))

def A_phi_MEG(A_phi_H2O, rho_H2O, epsilon_r_H2O, rho_MEG, epsilon_r_MEG):
    # Calcula Constante de Debye-Hückel para o MEG
    return ((rho_MEG/rho_H2O)**0.5 * (epsilon_r_H2O/epsilon_r_MEG)**1.5 * A_phi_H2O)

def A_phi_of_Mixing(rho_H2O_MEG,rho_H2O,epsilon_r_H2O_MEG,epsilon_r_H2O, T):
    #Calcula a constante de Debye-Hückel para mistura de solventes
    return ((rho_H2O_MEG/rho_H2O)**0.5)*((epsilon_r_H2O/epsilon_r_H2O_MEG)**1.5)*A_phi_H2O(T)
    
def pitzerParameters_MX_H2O(Beta_0, Beta_1, Cphi, T, z_M, Rm, Rx):
    #Calcula os parâmetros de Pitzer para a água
    Beta_0_MX_H2O = Beta_0 + (constants.K[0] * z_M**constants.Epsilon[0] + constants.K[1]*(Rm**constants.Epsilon[1]+Rx**constants.Epsilon[2])*z_M**constants.Epsilon[0]) * (T-constants.T_R)
    Beta_1_MX_H2O = Beta_1 + (constants.K[2] * z_M**constants.Epsilon[3] + constants.K[3]*(Rm**constants.Epsilon[4]+Rx**constants.Epsilon[5])*z_M**constants.Epsilon[3]) * (T-constants.T_R)
    Cphi_MX_H2O = Cphi + (constants.K[4] * z_M**constants.Epsilon[6] + constants.K[5]*(Rm**constants.Epsilon[7]+Rx**constants.Epsilon[8])*z_M**constants.Epsilon[6]) * (T-constants.T_R)    
    
    return Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O

def pitzerParameter_MX_MEG(nu_M, nu_X, A_phi_MEG, Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O, ln_gamma_MX_in_H2O, 
                           b_MX_in_H2O, b_MX_in_MEG, T, z_M, z_X, Rm, Rx, DeltaG_transf):
    #Calcula os parâmetros de Pitzer para o MEG
    Beta_0_MX_MEG = Beta_0_MX_H2O
    Cphi_MX_MEG = Cphi_MX_H2O
    
    def fobjetivoBeta1(beta1_MX):

        tamanho_v = len(beta1_MX)
        a = empty(tamanho_v)
        b = empty(tamanho_v)
        sol = 0.0
        
        for i in range(0,tamanho_v):
            a[i] = ln_gamma_MX(nu_M, nu_X, z_M, z_X, b_MX_in_MEG[i], A_phi_MEG[i], Beta_0_MX_MEG[i], beta1_MX[i], Cphi_MX_MEG[i]) 
            b[i] = ln_gamma_MX_in_MEG(nu_M, nu_X, ln_gamma_MX_in_H2O[i], b_MX_in_H2O[i], b_MX_in_MEG[i], T[i], DeltaG_transf)
            sol = sol + (a[i]-b[i])**2
        return sol
    
    def fobjetivoBeta1_Fsolve(beta1_MX):
        try:
            tamanho_v = len(beta1_MX)
            tamanho_t = len(T)
            a = empty(tamanho_v)
            b = empty(tamanho_v)
            if tamanho_t > 1:
                for i in range(0,tamanho_v):
                    a[i] = ln_gamma_MX(nu_M, nu_X, z_M, z_X, b_MX_in_MEG[i], A_phi_MEG[i], Beta_0_MX_MEG[i], beta1_MX[i], Cphi_MX_MEG[i]) 
                    b[i] = ln_gamma_MX_in_MEG(nu_M, nu_X, ln_gamma_MX_in_H2O[i], b_MX_in_H2O[i], b_MX_in_MEG[i], T[i], DeltaG_transf)
            if tamanho_t == 1:
                for i in range(0,tamanho_v):
                    a[i] = ln_gamma_MX(nu_M, nu_X, z_M, z_X, b_MX_in_MEG[i], A_phi_MEG[i], Beta_0_MX_MEG[i], beta1_MX[i], Cphi_MX_MEG[i]) 
                    b[i] = ln_gamma_MX_in_MEG(nu_M, nu_X, ln_gamma_MX_in_H2O[i], b_MX_in_H2O[i], b_MX_in_MEG[i], T, DeltaG_transf)
        except:
            a = ln_gamma_MX(nu_M, nu_X, z_M, z_X, b_MX_in_MEG, A_phi_MEG, Beta_0_MX_MEG, beta1_MX, Cphi_MX_MEG) 
            b = ln_gamma_MX_in_MEG(nu_M, nu_X, ln_gamma_MX_in_H2O, b_MX_in_H2O, b_MX_in_MEG, T, DeltaG_transf)
        return (a-b)**2

    tamanho_v = len(Beta_0_MX_MEG)
    lb = empty(tamanho_v)
    ub = empty(tamanho_v)
    for i in range(0, tamanho_v):
        lb[i] = -10.0
        ub[i] = 15.0
    
    Beta_1_MX_MEG_swarm, fopt = pso(fobjetivoBeta1, lb, ub, ieqcons=[], f_ieqcons=None, kwargs={}, swarmsize=100,
                 omega=0.5, phip=0.8, phig=0.8, maxiter=6000, minstep=1e-6, minfunc=1e-6, debug=False)
    
    #Beta_1_MX_MEG = minimize(fobjetivoBeta1,Beta_1_MX_MEG_swarm,method='Nelder-Mead', tol = 1e-6, options = {'xatol':1e-6, 'disp':True, 'maxiter': 10000})
    Beta_1_MX_MEG = fsolve(fobjetivoBeta1_Fsolve, Beta_1_MX_MEG_swarm)
    return Beta_0_MX_MEG, Beta_1_MX_MEG, Cphi_MX_MEG

def Pitzer_parameters_in_mixtures(Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O, Beta_1_MX_MEG, epsilon_r_H2O, epsilon_r_MEG, epsilon_r_H2O_MEG):
    # Calcula os parâmetros do Modelo de Pitzer em mistura H2O-MEG:
    # => beta0_MX_H2O_MEG, beta1_MX_H2O_MEG e C_phi_MX_H2O_MEG
    # => LORIMER, 1993
    Beta_0_MX_H2O_MEG = Beta_0_MX_H2O
    C_phi_MX_H2O_MEG = Cphi_MX_H2O
    tamanho_v = len(Beta_0_MX_H2O)
    A = empty(tamanho_v)
    B = empty(tamanho_v)
    Beta_1_MX_H2O_MEG = empty(tamanho_v)
    for i in range(0, tamanho_v):
        A[i] = log(Beta_1_MX_H2O[i]/Beta_1_MX_MEG[i]) / log(epsilon_r_H2O[i]/epsilon_r_MEG[i])
        B[i] = log(Beta_1_MX_H2O[i]) - A[i]*log(epsilon_r_H2O[i])
        Beta_1_MX_H2O_MEG[i] = exp(A[i]*log(epsilon_r_H2O_MEG[i]) + B[i])
    return Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, C_phi_MX_H2O_MEG

def Pitzer_parameters_in_mixtures_exc(Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O, Beta_1_MX_MEG, epsilon_r_H2O, epsilon_r_MEG, epsilon_r_H2O_MEG, x_MEG_SF):
    # Calcula os parâmetros do Modelo de Pitzer em mistura H2O-MEG:
    # => beta0_MX_H2O_MEG, beta1_MX_H2O_MEG e C_phi_MX_H2O_MEG
    # => LORIMER, 1993
    Beta_0_MX_H2O_MEG = Beta_0_MX_H2O
    C_phi_MX_H2O_MEG = Cphi_MX_H2O
    tamanho_v = len(x_MEG_SF)
    A = empty(tamanho_v)
    B = empty(tamanho_v)
    Beta_1_MX_H2O_MEG = empty(tamanho_v)
    A = log(Beta_1_MX_H2O/Beta_1_MX_MEG) / log(epsilon_r_H2O/epsilon_r_MEG)
    B = log(Beta_1_MX_H2O) - A*log(epsilon_r_H2O)
    for i in range(0, tamanho_v):
        Beta_1_MX_H2O_MEG[i] = exp(A*log(epsilon_r_H2O_MEG[i]) + B)
    return Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, C_phi_MX_H2O_MEG

def ln_gamma_MX(nu_M, nu_X, z_M, z_X, b_MX, A_phi, beta0_MX, beta1_MX, C_phi_MX):
    # Calcula o Coeficiente de Atividade pelo Modelo de Pitzer
    # => PITZER, 1973
    
    alpha = 2
    b_ = 1.2
    nu = nu_M + nu_X
    b_0 = 1 #Da tabela
    
    I_b = (abs(z_M*z_X) * b_MX * nu)/(2*b_0)
    sqrt_I_b = sqrt(I_b)
    try:
        tamanho_v = len(beta1_MX)
        if tamanho_v == 1:
            f_gamma_MX = -A_phi*((2/b_)*log(1 + b_*sqrt_I_b) + sqrt_I_b/(1 + b_*sqrt_I_b)) 
            B_gamma_MX = 2*(beta0_MX + beta1_MX*(1 - (1 + alpha*sqrt_I_b - 0.5*alpha**2*I_b)*exp(-alpha * sqrt_I_b))/(alpha**2*I_b))
            C_gamma_MX = 1.5*C_phi_MX
            ln_gamma_MX = (abs(z_M*z_X)*f_gamma_MX + 2*b_MX*nu_M*nu_X/nu*B_gamma_MX + (b_MX**2)*((2*(nu_M*nu_X)**1.5)*C_gamma_MX/nu))
        else:
            f_gamma_MX = empty(tamanho_v)
            B_gamma_MX = empty(tamanho_v)
            C_gamma_MX = empty(tamanho_v)
            ln_gamma_MX = empty(tamanho_v)
            for i in range(0, tamanho_v):
                f_gamma_MX[i] = -A_phi[i]*((2/b_)*log(1 + b_*sqrt_I_b[i]) + sqrt_I_b[i]/(1 + b_*sqrt_I_b[i])) 
                B_gamma_MX[i] = 2*(beta0_MX[i] + beta1_MX[i]*(1 - (1 + alpha*sqrt_I_b[i] - 0.5*alpha**2*I_b[i])*exp(-alpha * sqrt_I_b[i]))/(alpha**2*I_b[i]))
                C_gamma_MX[i] = 1.5*C_phi_MX[i]
                ln_gamma_MX[i] = (abs(z_M*z_X)*f_gamma_MX[i] + 2*b_MX[i]*nu_M*nu_X/nu*B_gamma_MX[i] + (b_MX[i]**2)*((2*(nu_M*nu_X)**1.5)*C_gamma_MX[i]/nu))
    except:
       f_gamma_MX = -A_phi*((2/b_)*log(1 + b_*sqrt_I_b) + sqrt_I_b/(1 + b_*sqrt_I_b)) 
       B_gamma_MX = 2*(beta0_MX + beta1_MX*(1 - (1 + alpha*sqrt_I_b - 0.5*alpha**2*I_b)*exp(-alpha * sqrt_I_b))/(alpha**2*I_b))
       C_gamma_MX = 1.5*C_phi_MX
       ln_gamma_MX = (abs(z_M*z_X)*f_gamma_MX + 2*b_MX*nu_M*nu_X/nu*B_gamma_MX + (b_MX**2)*((2*(nu_M*nu_X)**1.5)*C_gamma_MX/nu)) 
    
    return ln_gamma_MX

def ln_gamma_MX_in_MEG(nu_M, nu_X, ln_gamma_MX_in_H2O, b_MX_in_H2O, b_MX_in_MEG, T, DeltaG_transf):
    # Calcula o Coeficiente de Atividade do sal em MEG
    # => LORIMER, 1993
    nu = nu_M + nu_X
    return (ln_gamma_MX_in_H2O + log(b_MX_in_H2O/b_MX_in_MEG) - DeltaG_transf/(nu*constants.R * T))



def pitzerParameter_MX_MEG_SI(nu_M, nu_X, A_phi_MEG, Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O, ln_gamma_MX_in_H2O, 
                           b_MX_in_H2O, b_MX_in_MEG, T, z_M, z_X, DeltaG_transf):
    #Calcula os parâmetros de Pitzer para o MEG
    Beta_0_MX_MEG = Beta_0_MX_H2O
    Cphi_MX_MEG = Cphi_MX_H2O
    
    def fobjetivoBeta1(beta1_MX):

        a = ln_gamma_MX_SI(nu_M, nu_X, z_M, z_X, b_MX_in_MEG, A_phi_MEG, Beta_0_MX_MEG, beta1_MX, Cphi_MX_MEG) 
        b = ln_gamma_MX_in_MEG(nu_M, nu_X, ln_gamma_MX_in_H2O, b_MX_in_H2O, b_MX_in_MEG, T, DeltaG_transf)
        sol = (a-b)**2
        return sol
    
    def fobjetivoBeta1_Fsolve(beta1_MX):
        a = ln_gamma_MX_SI(nu_M, nu_X, z_M, z_X, b_MX_in_MEG, A_phi_MEG, Beta_0_MX_MEG, beta1_MX, Cphi_MX_MEG) 
        b = ln_gamma_MX_in_MEG(nu_M, nu_X, ln_gamma_MX_in_H2O, b_MX_in_H2O, b_MX_in_MEG, T, DeltaG_transf)
        return (a-b)**2

    lb = [-500.0]
    ub = [500.0]
    
    Beta_1_MX_MEG_swarm, fopt = pso(fobjetivoBeta1, lb, ub, ieqcons=[], f_ieqcons=None, kwargs={}, swarmsize=100,
                 omega=0.5, phip=0.8, phig=0.8, maxiter=6000, minstep=1e-6, minfunc=1e-6, debug=False)
   
    Beta_1_MX_MEG = minimize(fobjetivoBeta1,Beta_1_MX_MEG_swarm,method='Nelder-Mead', tol = 1e-6, options = {'xatol':1e-6, 'disp':True, 'maxiter': 10000})
    return Beta_0_MX_MEG, Beta_1_MX_MEG.x, Cphi_MX_MEG


def ln_gamma_MX_SI(nu_M, nu_X, z_M, z_X, b_MX, A_phi, beta0_MX, beta1_MX, C_phi_MX, len_Beta1_equals_to_1=True):
    # Calcula o Coeficiente de Atividade pelo Modelo de Pitzer
    # => PITZER, 1973
    
    alpha = 2
    b_ = 1.2
    nu = nu_M + nu_X
    b_0 = 1 #Da tabela
    C_gamma_MX = 1.5*C_phi_MX

    if len_Beta1_equals_to_1==True:
        I_b = (abs(z_M*z_X) * b_MX * nu)/(2*b_0)
        sqrt_I_b = sqrt(I_b)
        f_gamma_MX = -A_phi*((2/b_)*log(1 + b_*sqrt_I_b) + sqrt_I_b/(1 + b_*sqrt_I_b)) 
        B_gamma_MX = 2*(beta0_MX + beta1_MX*(1 - (1 + alpha*sqrt_I_b - 0.5*alpha**2*I_b)*exp(-alpha * sqrt_I_b))/(alpha**2*I_b))
        ln_gamma_MX = (abs(z_M*z_X)*f_gamma_MX + 2*b_MX*nu_M*nu_X/nu*B_gamma_MX + (b_MX**2)*((2*(nu_M*nu_X)**1.5)*C_gamma_MX/nu)) 
    else:
        tamanho_v = len(beta1_MX)
        B_gamma_MX = empty(tamanho_v)
        ln_gamma_MX = empty(tamanho_v)
        I_b = empty(tamanho_v)
        sqrt_I_b = empty(tamanho_v)
        f_gamma_MX = empty(tamanho_v)
        for i in range (0,tamanho_v):
            I_b[i] = (abs(z_M*z_X) * b_MX[i] * nu)/(2*b_0)
            sqrt_I_b[i] = sqrt(I_b[i])
            f_gamma_MX[i] = -A_phi[i]*((2/b_)*log(1 + b_*sqrt_I_b[i]) + sqrt_I_b[i]/(1 + b_*sqrt_I_b[i]))
            B_gamma_MX[i] = 2*(beta0_MX + beta1_MX[i]*(1 - (1 + alpha*sqrt_I_b[i] - 0.5*alpha**2*I_b[i])*exp(-alpha * sqrt_I_b[i]))/(alpha**2*I_b[i]))
            ln_gamma_MX[i] = (abs(z_M*z_X)*f_gamma_MX[i] + 2*b_MX[i]*nu_M*nu_X/nu*B_gamma_MX[i] + (b_MX[i]**2)*((2*(nu_M*nu_X)**1.5)*C_gamma_MX/nu))

    return ln_gamma_MX

def Pitzer_parameters_in_mixtures_SI(Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O, Beta_1_MX_MEG, epsilon_r_H2O, epsilon_r_MEG, epsilon_r_H2O_MEG):
    # Calcula os parâmetros do Modelo de Pitzer em mistura H2O-MEG:
    # => beta0_MX_H2O_MEG, beta1_MX_H2O_MEG e C_phi_MX_H2O_MEG
    # => LORIMER, 1993
    Beta_0_MX_H2O_MEG = Beta_0_MX_H2O
    C_phi_MX_H2O_MEG = Cphi_MX_H2O
    tamanho_v = len(epsilon_r_H2O_MEG)
    Beta_1_MX_H2O_MEG = empty(tamanho_v)
    A = log(Beta_1_MX_H2O/Beta_1_MX_MEG) / log(epsilon_r_H2O/epsilon_r_MEG)
    B = log(Beta_1_MX_H2O) - A*log(epsilon_r_H2O)
    for i in range(0, tamanho_v):
        Beta_1_MX_H2O_MEG[i] = exp(A*log(epsilon_r_H2O_MEG[i]) + B)
    return Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, C_phi_MX_H2O_MEG


def Pitzer_parameters_in_mixtures_exc_SI(Beta_0_MX_H2O, Beta_1_MX_H2O, Cphi_MX_H2O, Beta_1_MX_MEG, epsilon_r_H2O, epsilon_r_MEG, epsilon_r_H2O_MEG, x_MEG_SF):
    # Calcula os parâmetros do Modelo de Pitzer em mistura H2O-MEG:
    # => beta0_MX_H2O_MEG, beta1_MX_H2O_MEG e C_phi_MX_H2O_MEG
    # => LORIMER, 1993
    Beta_0_MX_H2O_MEG = Beta_0_MX_H2O
    C_phi_MX_H2O_MEG = Cphi_MX_H2O
    tamanho_v = len(x_MEG_SF)
    Beta_1_MX_H2O_MEG = empty(tamanho_v)
    A = log(Beta_1_MX_H2O/Beta_1_MX_MEG) / log(epsilon_r_H2O/epsilon_r_MEG)
    B = log(Beta_1_MX_H2O) - A*log(epsilon_r_H2O)
    for i in range(0, tamanho_v):
        Beta_1_MX_H2O_MEG[i] = exp(A*log(epsilon_r_H2O_MEG[i]) + B)
    return Beta_0_MX_H2O_MEG, Beta_1_MX_H2O_MEG, C_phi_MX_H2O_MEG