# Created by: Alejandra Risco (arisco@espol.edu.ec)

import FuncionesUNIFAC
import matplotlib.pyplot as plt
import numpy as np

print("********* UNIFAC Method *********\n")
l_comp = FuncionesUNIFAC.lcomp() #List of components on the database
l_structure,l_Antoine = FuncionesUNIFAC.l_infocomp()
l_op = list(range(len(l_comp)))

FuncionesUNIFAC.menu(l_comp)

n_comp1 = input("Component 1 (Most volatile):")
while FuncionesUNIFAC.val_ncomp(n_comp1,l_op):
    n_comp1 = input("\tCorrectly enter a # for Component 1:")
comp_1 = l_comp[int(n_comp1) - 1]

n_comp2 = input("Component 2 (Least volatile):")
while FuncionesUNIFAC.val_ncomp(n_comp2,l_op):
    n_comp2 = input("\tCorrectly enter a # for Component 1:")
comp_2 = l_comp[int(n_comp2) - 1]

sub, k_sub, Rk_sub, Qk_sub = FuncionesUNIFAC.subgroups()

#Parameters of each component
parametros_comp1 = FuncionesUNIFAC.parametersUNIFAC(l_structure,n_comp1,sub,k_sub,Rk_sub,Qk_sub)
parametros_comp2 = FuncionesUNIFAC.parametersUNIFAC(l_structure,n_comp2,sub,k_sub,Rk_sub,Qk_sub)

#r and q 
r1,q1 = FuncionesUNIFAC.ri_qi(parametros_comp1)
r2,q2 = FuncionesUNIFAC.ri_qi(parametros_comp2)

#e_k for subgroup
e_k1,k1 = FuncionesUNIFAC.e_ki(q1,parametros_comp1)
e_k2,k2 = FuncionesUNIFAC.e_ki(q2,parametros_comp2)

#t parameters Matrix
subgroups = FuncionesUNIFAC.unique_subgroups(parametros_comp1,parametros_comp2)
Temperatura = input("\nTemperature (C): ")
while not Temperatura.isdigit():
    Temperatura = input("Enter correctly a temperature in C:")
T = int(Temperatura) +273
print("\n")
M_t = FuncionesUNIFAC.subgroup_parameter(subgroups,T)

#Bi for each subgroup of each component
B_1 = FuncionesUNIFAC.Bi_comp(e_k1,k1,subgroups,M_t)
B_2 = FuncionesUNIFAC.Bi_comp(e_k2,k2,subgroups,M_t)

#Pxy Diagram
Psat_comp1 = FuncionesUNIFAC.P_sat(float(Temperatura),comp_1,n_comp1,l_Antoine)
Psat_comp2 = FuncionesUNIFAC.P_sat(float(Temperatura),comp_2,n_comp2,l_Antoine)


l_x,l_y,l_P = [],[],[]
for x1 in np.linspace(0,1,20): 
    if x1 == 0:
        y1 = 0
        P = float(Psat_comp2)
    elif x1 ==1:
        y1 = 1
        P = float(Psat_comp1)
    else:
        x2 = 1-x1
        coefy1,coefy2 = FuncionesUNIFAC.puntoUNIFAC(x1,x2,subgroups,r1,r2,q1,q2,B_1,B_2,k1,e_k1,k2,e_k2,M_t)
        P = x1*coefy1*float(Psat_comp1) + x2*coefy2*float(Psat_comp2)
        y1 = (x1*coefy1*float(Psat_comp1))/P
    l_x.append(x1)
    l_y.append(y1)
    l_P.append(P)

plt.plot(l_x,l_P,color='blue', label = 'Bubble line')
plt.plot(l_y,l_P,color='red', label = 'Dew line')
plt.xlabel('Composition (x1,y1) %s'%(comp_1.title()))
plt.ylabel('P [bar]')
s = 'Pxy Diagram %s/%s' %(comp_1.title(),comp_2.title())
plt.title(s)
plt.legend()
plt.show()