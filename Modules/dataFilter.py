from numpy import empty, mean, nan

def filter(x, y, z, T_All, T_filter):
    T_All_size = len(T_All)
    T_filter_size = len(T_filter)
    x_L, y_L, z_L = [], [], []
    x_filter, y_filter, z_filter = empty(T_filter_size), empty(T_filter_size), empty(T_filter_size)
    cont = 0
    for i in range(1, T_All_size):
    
        if T_All[i-1] == T_All[i]:
            x_L.append(x[i-1])
            y_L.append(y[i-1])
            z_L.append(z[i-1])
        
        else:
            try:
                if x_L[0] != nan:
                    x_filter[cont] = mean(x_L) 
                    y_filter[cont] = mean(y_L) 
                    z_filter[cont] = mean(z_L) 

                    cont+=1
                
                    x_L.clear()
                    y_L.clear()
                    z_L.clear()
                else:
                    x_filter[cont] = x[i-1] 
                    y_filter[cont] = y[i-1]
                    z_filter[cont] = z[i-1] 
                    cont+=1  
            except:
                x_filter[cont] = x[i-1] 
                y_filter[cont] = y[i-1]
                z_filter[cont] = z[i-1] 
                cont+=1

        if T_All_size == i+1:
        
            x_L.append(x[i])
            y_L.append(y[i])
            z_L.append(z[i])

            x_filter[cont] = mean(x_L)
            y_filter[cont] = mean(y_L)
            z_filter[cont] = mean(z_L)

    return x_filter, y_filter, z_filter

def filter_1(x, T_All_treated, T_filter):
    list_all = [[] for _ in range(0,len(T_filter))]
    cont = 0
    for i in range(1, len(T_All_treated)):
        if T_All_treated[i-1] == T_All_treated[i]:
            list_all[cont].append(x[i-1])
        else:
            list_all[cont].append(x[i-1])
            cont+=1
        if len(T_All_treated) == i+1:
            list_all[cont].append(x[i])
    return list_all