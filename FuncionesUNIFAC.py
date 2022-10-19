# Created by: Alejandra Risco (arisco@espol.edu.ec)

def lcomp():
    archivo = open("Properties_Comp.csv", "r")
    archivo.readline()
    lcomp = []
    for linea in archivo:
        l= linea.strip().split(";")
        lcomp.append(l[0])
    archivo.close()
    return lcomp

def l_infocomp():
    archivo = open("Properties_Comp.csv", "r")
    archivo.readline()
    lstructure = []
    lAntoine = []
    for linea in archivo:
        l= linea.strip().split(";")
        lAntoine.append(tuple(l[1:-1]))
        lstructure.append(l[-1])
    archivo.close()
    return lstructure,lAntoine

def menu(l_comp):
    print("Available components:")
    rango = list(range(0, len(l_comp), 4))
    for i in rango:
        print("{:>3}. {:20} {:>3}. {:20} {:>3}. {:20} {:>3}. {:20}".format(i+1, l_comp[i].title(),i+2,l_comp[i+1].title(),i+3,l_comp[i+2].title(),i+4,l_comp[i+3].title()))
    print("\n")
    
def val_ncomp(ncomp,l_op):
    if not ncomp.isdigit():
        return True
    elif ncomp.isdigit():
        if int(ncomp) - 1 not in l_op:
            return True
        else:
            return False

def coef_Antoine(n_comp,l_Antoine):
    A = float(l_Antoine[int(n_comp)-1][0].replace(",","."))
    B = float(l_Antoine[int(n_comp)-1][1].replace(",","."))
    C = float(l_Antoine[int(n_comp)-1][2].replace(",","."))
    T_min = int(l_Antoine[int(n_comp)-1][3].split(",")[0])
    T_max = int(l_Antoine[int(n_comp)-1][4].split(",")[0])
    return [A,B,C,list(range(T_min,T_max+1))]

def val_numero(x):
    while x.isalpha():
        x = input("\tEnter correctly a number:")
    if "," in x:
        x=x.replace(",",".")
    return float(x)

def P_sat(T,comp,n_comp,l_Antoine):
    A,B,C,rangeT = coef_Antoine(n_comp,l_Antoine)
    T1 = round(T,0)
    if T1 in rangeT:
        Psat_comp = 10**(A-((B)/(T+C))) #mmHg
        Psat_comp = Psat_comp*(1.01325/760) #bar
    else:
        Psat_comp = input("Antoine coefficients for %s are not valid for the T of the system\nEnter its Psat [bar] a %.2f C: "%(comp,T)) #bar
        Psat_comp = val_numero(Psat_comp)
    return Psat_comp


def subgroups():
    archivo = open("Ksubgroup_Parameters.csv","r")
    archivo.readline()
    subgrupos = []
    k_sub = []
    Rk_sub = []
    Qk_sub = []
    for linea in archivo:
        data = linea.split(",")
        subgrupo = data[1]
        subgrupos.append(subgrupo)
        k_sub.append(int(data[2]))
        Rk_sub.append(float(data[3]))
        Qk_sub.append(float(data[4]))
    return subgrupos,k_sub,Rk_sub,Qk_sub

def parametersUNIFAC(l_s,n_comp,sub_l,k_l,Rk_l,Qk_l):
    comp = l_s[int(n_comp)-1].split('-')
    parametros_comp = []
    for subgrupo in comp:
        if "|" in subgrupo:
            vk,sub = int(subgrupo.split("|")[0]),subgrupo.split("|")[1]
        else:
            vk,sub = 1,subgrupo
        ind = sub_l.index(sub)
        k_subgrupo = k_l[ind]
        Rk_subgrupo = Rk_l[ind]
        Qk_subgrupo =Qk_l[ind]
        tupla = (sub,vk,k_subgrupo,Rk_subgrupo,Qk_subgrupo)
        parametros_comp.append(tupla)
    return parametros_comp

def ri_qi(parametros_comp):
    ri,qi = 0,0
    for subi in parametros_comp:
        vki = subi[1]
        Rki = subi[3]
        Qki = subi[4]
        ri += vki*Rki
        qi += vki*Qki
    return ri,qi

def e_ki(qi,parametros_comp):
    e_ki = []
    k_i = []
    for subi in parametros_comp:
        ki = subi[0]
        vki = subi[1]
        Qki = subi[4]
        eki = (vki*Qki)/qi
        e_ki.append(eki)
        k_i.append(ki)
    return e_ki,k_i

def unique_subgroups(parametros1,parametros2):
    subgrupos = []
    for i in range(len(parametros1)):
        sub = parametros1[i][0]
        if sub not in subgrupos:
            subgrupos.append(sub)
    for i in range(len(parametros2)):
        sub = parametros2[i][0]
        if sub not in subgrupos:
            subgrupos.append(sub)
    subgrupos = list(set(subgrupos))
    return subgrupos

def group_subgroups(subgrupo):
    archivo = open("Ksubgroup_Parameters.csv", "r")
    archivo.readline()
    grupos_principales = []
    sub_gruposprincipales = []
    for linea in archivo:
        data = linea.split(",")
        grupo = data[0]
        if grupo not in grupos_principales:
            grupos_principales.append(grupo)
            sub_gruposprincipales.append([data[1]])
        else:
            ind = grupos_principales.index(grupo)
            sub_gruposprincipales[ind].append(data[1])
    for i in sub_gruposprincipales:
        if subgrupo in i:
            ind = sub_gruposprincipales.index(i)
            grupo = grupos_principales[ind]
            return grupo

def interaction_parameter(grupo1,grupo2):
    archivo = open("Interaction_Parameters.csv","r",encoding='utf-8')
    l_grupo2 = archivo.readline().strip().split(",")
    ind = l_grupo2.index(grupo2)
    for linea in archivo:
        data =linea.split(",")
        if grupo1 == data[0]:
            parametro_inter = data[ind]
    return parametro_inter

def subgroup_parameter(subgrupos,Temperatura_K):
    import numpy as np
    import math
    Matriz_parametrost = np.zeros((len(subgrupos),len(subgrupos)),dtype=float)
    for fila in range(len(subgrupos)):
        for columna in range(len(subgrupos)):
            grupo_columna = group_subgroups(subgrupos[columna])
            grupo_fila = group_subgroups(subgrupos[fila])
            parametro_a = float(interaction_parameter(grupo_fila,grupo_columna))
            Matriz_parametrost[fila,columna] = math.exp(-parametro_a/Temperatura_K)
    return Matriz_parametrost

def Bi_comp(e_ki,k_i,subgrupos,Matriz_parametrost):
    Bi_l = []
    for k in subgrupos:
        ind = subgrupos.index(k)
        Matriz_parametrosti = Matriz_parametrost[:,ind]
        Bi = 0
        for i in range(len(subgrupos)):
            if subgrupos[i] in k_i:
                ind = k_i.index(subgrupos[i])
                ei = e_ki[ind]
            else:
                ei = 0
            idx = subgrupos.index(subgrupos[i])
            ti = Matriz_parametrosti[idx]
            Bi += ei * ti
        Bi_l.append(Bi)
    return Bi_l

def tetha_subgroup(subgrupos,q1,q2,ki1,e_ki1,ki2,e_ki2,x1,x2):
    tetha_l = []
    for k in subgrupos:
        if k in ki1:
            idx = ki1.index(k)
            ek1 = e_ki1[idx]
        else:
            ek1 = 0
        if k in ki2:
            idx = ki2.index(k)
            ek2 = e_ki2[idx]
        else:
            ek2 = 0
        tetha_k = (x1*q1*ek1+x2*q2*ek2)/(x1*q1+x2*q2)
        tetha_l.append(tetha_k)
    return tetha_l

def S_subgroup(subgrupos,tetha_l,Matriz_parametrost):
    import numpy as np
    S_l = []
    for k in subgrupos:
        ind = subgrupos.index(k)
        Matriz_parametrost_k = Matriz_parametrost[:,ind]
        tetha_k = np.array(tetha_l)
        S_k = np.dot(Matriz_parametrost_k,tetha_k)
        S_l.append(S_k)
    return S_l

def J_L_comp(r1,r2,x1,x2,q1,q2):
    J1 = (r1)/(r1*x1+r2*x2)
    L1 = (q1)/(q1*x1+q2*x2)
    J2 = (r2) / (r1 * x1 + r2 * x2)
    L2 = (q2) / (q1 * x1 + q2 * x2)
    return (J1,L1),(J2,L2)

def ln_yi_C(Ji,Li,qi):
    import numpy as np
    ln_yi = 1 - Ji + np.log(Ji) - 5*qi*(1-(Ji/Li)+np.log((Ji/Li)))
    return ln_yi

def ln_yi_R(qi,Bi,Si,tethai,subgrupos,ki,eki):
    import numpy as np
    suma = 0
    for k in subgrupos:
        if k in ki:
            idx = ki.index(k)
            ei_k = eki[idx]
        else:
            ei_k = 0
        ind = subgrupos.index(k)
        Bi_k = Bi[ind]
        Si_k = Si[ind]
        tethai_k = tethai[ind]
        ec = tethai_k*(Bi_k/Si_k)-ei_k*np.log((Bi_k/Si_k))
        suma += ec
    ln_yi = qi*(1-suma)
    return ln_yi

def yi(ln_yi_C,ln_yi_R):
    import math
    lnyi_total = ln_yi_C + ln_yi_R
    yi = math.exp(lnyi_total)
    return yi


def puntoUNIFAC(x1,x2,subgrupos,r1,r2,q1,q2,B_1,B_2,k1,e_k1,k2,e_k2,M_t):
  tetha_sub = tetha_subgroup(subgrupos,q1,q2,k1,e_k1,k2,e_k2,x1,x2)
  S_sub = S_subgroup(subgrupos,tetha_sub,M_t)
  J_L_comp1,J_L_comp2 = J_L_comp(r1,r2,x1,x2,q1,q2)
  J1,L1 = J_L_comp1[0],J_L_comp1[1]
  J2,L2 = J_L_comp2[0],J_L_comp2[1]
  
  #Component 1
  ln_y1_C =ln_yi_C(J1,L1,q1)
  ln_y1_R = ln_yi_R(q1,B_1,S_sub,tetha_sub,subgrupos,k1,e_k1)
  
  #Component 2
  ln_y2_C = ln_yi_C(J2,L2,q2)
  ln_y2_R = ln_yi_R(q2,B_2,S_sub,tetha_sub,subgrupos,k2,e_k2)
  
  #Activity Coefficient
  y1 = yi(ln_y1_C,ln_y1_R)
  y2 = yi(ln_y2_C,ln_y2_R)
  return y1,y2