from numpy import amax, sum, mean, empty
# R² calculation
def coefficient_of_determination(y, ycalc):
    mean = sum(y)/len(y)
    residue = (ycalc - mean)**2
    SQ_residue = sum(residue)
    variance = (y - mean)**2
    SQ_Total = sum(variance)
    R_2 = SQ_residue/SQ_Total
    return R_2

# R² calculation
def coefficient_of_determination2(y, ycalc):
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
    # Mean absolute deviation
    AAD = mean(AD)
    # Maximum absolute deviation
    maxAD = amax(AD)
    # Mean Absolute Relative Deviation
    AARD = mean(ARD)
    # Maximum Absolute Relative Deviation
    maxARD = amax(ARD)
    return AAD, maxAD, AARD, maxARD


def coefficient_of_determination_JAFO(y, ycalc):
    mean = sum(y)/len(y)
    # Total Sum of Squares (TSS)
    TSS = sum((y-mean)**2)
    # Residual Sum of Squares (RSS)
    RSS = sum((y-ycalc)**2)
    R_2 = (TSS - RSS)/TSS
    return R_2   
