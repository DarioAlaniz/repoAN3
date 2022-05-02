# https://scikit-rf.readthedocs.io/en/latest/ libreria para leer parametros S

import numpy as np
import skrf as  rf 
import matplotlib.pylab as plt
from sympy import deg

def polar(number):
    return (np.abs(number),np.angle(number,deg=True))

z0=50       #[ohm]
f=1.4e9     #Cambiar a 1600[Gz] Grupo 14!! 1.4 es de ejemplo para ver si anda bien el programa
Vce = 1.0   #1 a 4 [V] con paso de 0.5
Ic  = 42    #18-22-28-35-37-42-48-50 [mA]

########################################
# Obtencion de parametros S
########################################
path = "SPAR\BFP640\BFP640_VCE_{}V_IC_{}mA.s2p".format(Vce,Ic)
s2p = rf.Network(path)
ind = np.where(s2p.f == f)[0]   #Encuentro el indice donde se encuentra la freq elegida para obtener los parametros S
param = s2p.s[ind][0]           #Parametros S encontrados expresados en Real Imaginario 
print(print("Parametros S obtenidos a f={}[GHz] y polarizacion Ic={}[mA],Vce={}[V]".format(format(f,'.1E'),Ic,Vce).center(100,"_")))
s11 = param[0][0]
s12 = param[0][1]
s21 = param[1][0]
s22 = param[1][1]
for i in range(len(param)):
    for j in range(len(param)):
        print("S_{}{}={}".format(i+1,j+1,polar(param[i][j])))
########################################
# Analisis de estabilidad
########################################
# k = s2p.stability[ind]
print('Estabilidad'.center(100,"_"))
D = np.abs((s11*s22)-(s12*s21))
k = (1-(np.abs(s22)**2)-(np.abs(s11)**2) + (D**2) ) / (2*np.abs(s12*s21))
print('Factor de Rollet (k)={}'.format(k))
print('Coeficiente delta: {}'.format(D))
if(k<1):
    print('Condicionalmente estable, elegir otra polarizacion')
    exit()
else: 
    print("Incondicionalmente Estable, se puede seguir con el diseño para maxima ganancia.")
########################################
# Calculos auxiliares
########################################
B1 = 1 + (np.abs(s11)**2) - np.abs(s22)**2 - (D**2)
C1 = s11 - (D*np.conjugate(s22)) 
B2 = 1 + (np.abs(s22)**2) - np.abs(s11)**2 - (D**2)
C2 = s22 - (D*np.conjugate(s11))
print(B1,polar(C1),B2,polar(C2))
########################################
# Calculos de Coeficientes de reflexion 
########################################
print('Coeficientes de Relexion'.center(100,'_'))
print('Entrada: ')
r_max_in =  ((B1 - np.sqrt((B1**2)-(4*(np.abs(C1)**2))))/(2*(np.abs(C1)**2))) * np.conjugate(C1)
print(polar(r_max_in))
print(r_max_in)
print('Salida: ')
r_max_out = ((B2 - np.sqrt((B2**2)-(4*(np.abs(C2)**2))))/(2*(np.abs(C2)**2)))*np.conjugate(C2)
print(polar(r_max_out))
print(r_max_out)
########################################
# Calculos de Impedancias de entrada y salida
########################################
print('Impedancias de entrada y salida'.center(100,'_'))
print('Entrada: ')
z_in_serie = z0 * ((1+r_max_in)/(1-r_max_in))
print(z_in_serie)
print('Salida: ')
z_out_serie = z0 * ((1+r_max_out)/(1-r_max_out))
print(z_out_serie)