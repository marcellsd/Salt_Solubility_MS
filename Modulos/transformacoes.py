# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 17:09:45 2021

@author: marcell
"""
from Modulos import constants

def fracao_massica_para_molalidade(MM_MX, w_MX):
    return (w_MX/MM_MX)/(1-w_MX)


def x_MEG_SF(w_H2O_SF):
    return (((1-w_H2O_SF)/constants.MM_MEG)/((w_H2O_SF/constants.MM_H2O)+((1-w_H2O_SF)/constants.MM_MEG)))

def w_MEG_SF(x_MEG_SF):
    return (x_MEG_SF*constants.MM_MEG)/(x_MEG_SF*constants.MM_MEG + (1-x_MEG_SF)*constants.MM_H2O)
