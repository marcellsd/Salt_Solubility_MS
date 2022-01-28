from numpy import empty, mean, nan

def filtro(x, y, z, T_All, T_filtro):
    T_All_size = len(T_All)
    T_filtro_size = len(T_filtro)
    x_L, y_L, z_L = [], [], []
    x_filtro, y_filtro, z_filtro = empty(T_filtro_size), empty(T_filtro_size), empty(T_filtro_size)
    cont = 0
    for i in range(1, T_All_size):
    
        if T_All[i-1] == T_All[i]:
            x_L.append(x[i-1])
            y_L.append(y[i-1])
            z_L.append(z[i-1])
        
        else:
            try:
                if x_L[0] != nan:
                    x_filtro[cont] = mean(x_L) 
                    y_filtro[cont] = mean(y_L) 
                    z_filtro[cont] = mean(z_L) 

                    cont+=1
                
                    x_L.clear()
                    y_L.clear()
                    z_L.clear()
                else:
                    x_filtro[cont] = x[i-1] 
                    y_filtro[cont] = y[i-1]
                    z_filtro[cont] = z[i-1] 
                    cont+=1  
            except:
                x_filtro[cont] = x[i-1] 
                y_filtro[cont] = y[i-1]
                z_filtro[cont] = z[i-1] 
                cont+=1

        if T_All_size == i+1:
        
            x_L.append(x[i])
            y_L.append(y[i])
            z_L.append(z[i])

            x_filtro[cont] = mean(x_L)
            y_filtro[cont] = mean(y_L)
            z_filtro[cont] = mean(z_L)

    return x_filtro, y_filtro, z_filtro

def filtro_1(x, T_All_tratado, T_filtro):
    list_all = [[] for _ in range(0,len(T_filtro))]
    cont = 0
    for i in range(1, len(T_All_tratado)):
        if T_All_tratado[i-1] == T_All_tratado[i]:
            list_all[cont].append(x[i-1])
        else:
            list_all[cont].append(x[i-1])
            cont+=1
        if len(T_All_tratado) == i+1:
            list_all[cont].append(x[i])
    return list_all