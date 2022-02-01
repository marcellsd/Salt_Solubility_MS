from numpy import amax, sum, mean, empty
#Calculo R2
def coef_determinacao(y, ycalc):
    mean = sum(y)/len(y)
    residuo = (ycalc - mean)**2
    SQ_residuo = sum(residuo)
    variancia = (y - mean)**2
    SQ_Total = sum(variancia)
    R_2 = SQ_residuo/SQ_Total
    return R_2


def coef_determinacao2(y, ycalc):
    SSres_v = empty(len(y))
    for i in range(len(SSres_v)):
        SSres_v[i] = (y[i]-ycalc[i])**2
    SSres = sum(SSres_v)
    mean = sum(y)/len(y)
    Total = (y - mean)**2
    SStotal = sum(Total)
    R_2 = 1 - SSres/SStotal
    return R_2

def DAM_DMR(yexp, ycalc):
    
    AD = empty(len(yexp))
    ARD = empty(len(yexp))
    for i in range(len(AD)):
        AD[i] = abs(yexp[i]-ycalc[i])
        ARD[i] = abs((yexp[i]-ycalc[i])/yexp[i])
    # Desvio Absoluto Médio
    AAD = mean(AD)
    # Desvio Absoluto Máximo
    maxAD = amax(AD)
    # Desvio Relativo Absoluto Médio
    AARD = mean(ARD)
    # Desvio Relativo Absoluto Máximo
    maxARD = amax(ARD)
    return AAD, maxAD, AARD, maxARD


def coef_determinacao_JAFO(y, ycalc):
    mean = sum(y)/len(y)
    # Soma de Quadrados Total (TSS)
    TSS = sum((y-mean)**2)
    # Soma de Quadrados dos Resíduos (RSS)
    RSS = sum((y-ycalc)**2)
    R_2 = (TSS - RSS)/TSS
    return R_2   
