# https://scikit-rf.readthedocs.io/en/latest/ libreria para leer parametros S
import numpy as np
import skrf as  rf 
import matplotlib.pylab as plt
import sympy as sym     #simbolica
sym.init_printing()
from microtiraDesing import *
import cmath

def cart2polar(number,deg=True):
    return (np.abs(number),np.angle(number,deg=deg))

def polar2cart(mod,phi):
    return  complex(mod * np.cos(np.deg2rad(phi)) , mod * np.sin(np.deg2rad(phi)))

def serie2parallel(z):
    a = np.real(z) * (1 + (np.imag(z)/np.real(z))**2)
    b = np.imag(z) * (1 + (np.real(z)/np.imag(z))**2)
    return a,b

def adapt(z):
    r_p,x_p = serie2parallel(z)
    if(x_p>0): #Inductor
        print('No se necesita usar el transformador lamda/4 para mover la reactancia ya que es inductiva')
        z_lamba4 = np.sqrt(z0*r_p)
    else:
        print('Se necesita usar el transformador lamda/4 para mover la reactancia ya que es capacitiva')
        z_lamba4 = np.sqrt(z0*np.real(z)) #se trabaja directamente con el modelo serie ya al usar el lamda/4 se tiene la inductancia en paralelo
        x_p    = (z_lamba4**2) / np.imag(z)       #pasa de ser capacitiva a inductiva 
    C  = 1/(2*np.pi* f * np.abs(x_p)) 
    return z_lamba4,C,r_p,x_p

def length_C(z0,f,er_p,C):
    lamda   = 300/(f/1e6)
    lamda_p = lamda/np.sqrt(er_p)
    beta    = (2*np.pi)/(lamda_p)
    X       = 1/(2*np.pi*f*C)
    d       = np.arctan(z0/X)/beta
    return d
def length_lamda4(f,er_p):
    lamda       =  300/(f/1e6)
    lamda_p     =  lamda/np.sqrt(er_p)
    return  (lamda_p/4)
##########################################################################
# Para calcular los circuilos de estabilidad
sqabs = lambda x: np.square(np.absolute(x))
def rad_center(s11,s12,s21,s22,delta):
    # source
    rs = np.absolute((s12 * s21)/(sqabs(s11) - sqabs(delta)))
    cs = np.conj(s11 - delta*np.conj(s22))/(sqabs(s11) - sqabs(delta))
    # load
    rl = np.absolute((s12 * s21)/(sqabs(s22) - sqabs(delta)))
    cl = np.conj(s22 - delta*np.conj(s11))/(sqabs(s22) - sqabs(delta))
    return rs,cs,rl,cl

def calc_circle(c, r):
    theta = np.linspace(0, 2*np.pi, 1000)
    return c + r*np.exp(1.0j*theta)
###########################################################################

Z_out , Z_in    = sym.symbols('Z_out,Z_in')
r_in , r_out    = sym.symbols('\sigma_in,\sigma_out')
MGS , MAG       = sym.symbols('MGS_db,MAG_db')

############################################################################
# Parametros de Diseño
z0  = 50        #[ohm]
zL  = 50
f   = 1.6e9     #Cambiar a 1600[Gz] Grupo 14!! 1.4 es de ejemplo para ver si anda bien el programa, 1.8 de test
Vce = 1.0       #1 a 4 [V] con paso de 0.5 ,1 v de test
Ic  = 37        #18-22-28-32-37-42-48-50 [mA] , 60 mA de test
#Mejor polarizacion obtenida 1.0v y 37mA. Dan valores de W entre 1 y 10
########################################
# Obtencion de parametros S
########################################
# path = "SPAR\BFP450\BFP450_w_noise_VCE_{}V_IC_{}mA.s2p".format(Vce,Ic) #path de test
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
        print("S_{}{}={}".format(i+1,j+1,cart2polar(param[i][j])))
########################################
# Analisis de estabilidad
########################################
# k = s2p.stability[ind]
print('Estabilidad'.center(100,"_"))
D = (s11*s22)-(s12*s21)
k = (1-(np.abs(s22)**2)-(np.abs(s11)**2) + (np.abs(D)**2) ) / (2*np.abs(s12*s21))
print('Factor de Rollet (k)={}'.format(k))
print('Coeficiente delta: {}'.format(D),'mod => {}'.format(np.abs(D)))
if(k<1):
    print('Condicionalmente estable, verificar circulos de estabilidad')
else: 
    print("Incondicionalmente Estable, se puede seguir con el diseño para maxima ganancia.")
# Plot circulos de estabilidad
# Circulos de estabilidad
# rs,cs,rl,cl = rad_center(s11,s12,s21,s22,D)
# ns = rf.Network(name='f='+str(f/1.e+9)+', S1', s=calc_circle(cs, rs))
# ns.plot_s_smith(draw_labels=True)
# nl = rf.Network(name='f='+str(f/1.e+9)+', S2', s=calc_circle(cl, rl))
# nl.plot_s_smith(draw_labels=True)
# plt.show()
# a = input('1: Para seguir, 0: Para salir ')
# if(a=='0'):
#     exit()
# s2p.plot_s_db()
# plt.show()
########################################
# Calculos auxiliares
########################################
B1 = 1 + (sqabs(s11)) - sqabs(s22) - (sqabs(D))
C1 = s11 - (D*np.conjugate(s22)) 
B2 = 1 + (sqabs(s22)) - sqabs(s11) -  sqabs(D)
C2 = s22 - (D*np.conjugate(s11))
print('B1={}, B2={}, C1={}, C2={} '.format(B1,B2,C1,C2))
########################################
# Calculos de Coeficientes de reflexion 
########################################
print('Coeficientes de Relexion'.center(100,'_'))
sym.print_latex(r_in)
if (B1 > 0):
    r_max_in =  np.abs((B1 - cmath.sqrt((B1**2)-(4*(sqabs(C1)))))/(2*(np.abs(C1)))) 
else :
    r_max_in =  np.abs((B1 + cmath.sqrt((B1**2)-(4*(sqabs(C1)))))/(2*(np.abs(C1))))

r_max_in = polar2cart(r_max_in,np.angle(C1,deg=True))
print(r_max_in)
print(cart2polar(r_max_in))

sym.print_latex(r_out)
if(B2>0):
    r_max_out = np.abs((B2 - cmath.sqrt((B2**2)-(4*(np.abs(C2)**2))))/(2*(np.abs(C2))))
else:
    r_max_out = np.abs((B2 + cmath.sqrt((B2**2)-(4*(np.abs(C2)**2))))/(2*(np.abs(C2))))

r_max_out = polar2cart(r_max_out,np.angle(C2,deg=True))
print(r_max_out)
print(cart2polar(r_max_out))

########################################
# Calculos de Impedancias de entrada y salida
########################################
print('Impedancias de entrada y salida'.center(100,'_'))
sym.print_latex(Z_in)
z_in_serie = z0 * ((1+r_max_in)/(1-r_max_in))
print(z_in_serie)
sym.print_latex(Z_out)
z_out_serie = z0 * ((1+r_max_out)/(1-r_max_out))
print(z_out_serie)
########################################
# Calculos de Ganancia
########################################
print('Ganancias'.center(100,'_'))
# Maximun stable Gain
mgs_db = 10 *np.log10(np.abs(s21)/np.abs(s12))
sym.print_latex(MGS)
print(mgs_db)
# Maximun aviable Gain
mag_db = 10 *np.log10((np.abs(s21)/np.abs(s12)) * (k - np.sqrt((k**2)-1)) )
sym.print_latex(MAG)
print(mag_db)
########################################
# Adaptacion en la entrada
########################################
print('Adaptacion en la entrada'.center(100,'_'))
z_lamba4_in,C_in,R_in_p , X_in_p = adapt(z_in_serie)
# print(z_lamba4_in,C_in,R_in_p , X_in_p)
# Modelo paralelo de entrada
if(np.angle(X_in_p,deg=True)>=0 and np.angle(X_in_p,deg=True)<180): # si es mayor a cero es que trabajo con el modelo paralelo y no serie
    print('Modelo paralelo de entrada'.center(100,'_'))
    Rinp, Xinp = sym.symbols('R_in(p) , X_in(p)')
    R_in_p , X_in_p = serie2parallel(z_in_serie) 
    sym.print_latex(Rinp)
    print(R_in_p)
    sym.print_latex(Xinp)
    print(X_in_p)
else:
    print('reactancia desplazada lamda/4 es: {}'.format(X_in_p))

print('Transformador lamda/4 de entrada'.center(100,'_'))
zlamba4in   = sym.symbols('Z_\lamda/4(in)')
sym.print_latex(zlamba4in)
print(z_lamba4_in)
print('Capacitor de Entrada'.center(100,'_'))
print('Se necesita sintetizar un capacitor de {}[F]'.format(C_in))
########################################
# Adaptacion en la salida
########################################
print('Adaptacion en la salida'.center(100,'_'))
z_lamba4_out,C_out,R_out_p, X_out_p  = adapt(z_out_serie)
# print(z_lamba4_out,C_out,R_out_p, X_out_p,np.angle(X_out_p,deg=True))
# Modelo paralelo de salida
if(np.angle(X_out_p,deg=True)>=0 and np.angle(X_out_p,deg=True)<180):# si es mayor a cero es que se trabajo con el modelo paralelo y no serie
    print('Modelo paralelo de salida'.center(100,'_'))
    Routp, Xoutp = sym.symbols('R_out(p) , X_out(p)')
    sym.print_latex(Routp)
    print(R_out_p)
    sym.print_latex(Xoutp)
    print(X_out_p)
else:
    print('reactancia desplazada lamda/4 es: {}'.format(X_out_p))

print('Transformador lamda/4 de salida'.center(100,'_'))
zlamba4out   = sym.symbols('Z_\lamda/4(out)')
sym.print_latex(zlamba4out)
print(z_lamba4_out)
print('Capacitor de salida'.center(100,'_'))
print('Se necesita sintetizar un capacitor de {}[F]'.format(C_out))

########################################
# Microtiras 
########################################
er= 4.5
H = 1.66  #mm datos de ejemplo
t = 0.04  #mm datos de ejemplo
tipo = 0  # 1 Hammerstad, 0 Wheeler

# Transformador lamba/4 de entrada
z = z_lamba4_in   #[ohm] impedancia a diseñar
print('Diseño de microtira para Transformador lamba/4 de entrada con z0= {}'.format(z).center(100,'_'))
A,B,W,ratio,er_p,z0_p   = microtira(tipo,t,H,er,z)
d_lamba4_in             =  length_lamda4(f,er_p)
print('Longitud de microtira: {}[m]'.format(d_lamba4_in))

# Transformador lamba/4 de salida
z = z_lamba4_out   #[ohm] impedancia a diseñar
print('Diseño de microtira para Transformador lamba/4 de salida con z0= {}'.format(z).center(100,'_'))
A,B,W,ratio,er_p,z0_p   = microtira(tipo,t,H,er,z)
d_lamba4_out            =  length_lamda4(f,er_p)
print('Longitud de microtira: {}[m]'.format(d_lamba4_out))

# Capacitor de entrada
z = 50   #[ohm] impedancia a diseñar Consultar si usamos z_lamba4_in o 50 ohm
print('Diseño de microtira para capacitor de entrada con z0= {}'.format(z).center(100,'_'))
A,B,W,ratio,er_p,z0_p   = microtira(tipo,t,H,er,z)
d_Cin = length_C(z0_p,f,er_p,C_in)
print('Longitud de microtira: {}[m]'.format(d_Cin))


# Capacitor de salida
z = 50   #[ohm] impedancia a diseñar
print('Diseño de microtira para capacitor de salida con z0= {}'.format(z).center(100,'_'))
A,B,W,ratio,er_p,z0_p   = microtira(tipo,t,H,er,z)
d_Cout = length_C(z0_p,f,er_p,C_out)
print('Longitud de microtira: {}[m]'.format(d_Cout))

# Transformador lamba/4 de desacople 1
z = 25   #[ohm] impedancia a diseñar
print('Diseño de microtira para capacitor de desacople con z0= {}'.format(z).center(100,'_'))
A,B,W,ratio,er_p,z0_p   = microtira(tipo,t,H,er,z)
d_lamba4_descople_CA             =  length_lamda4(f,er_p)
print('Longitud de microtira: {}[m]'.format(d_lamba4_descople_CA))

# Transformador lamba/4 de desacople 2
z = 95   #[ohm] impedancia a diseñar
print('Diseño de microtira para Inductor choque con z0= {}'.format(z).center(100,'_'))
A,B,W,ratio,er_p,z0_p   = microtira(tipo,t,H,er,z)
d_lamba4_descople_CC             =  length_lamda4(f,er_p)
print('Longitud de microtira: {}[m]'.format(d_lamba4_descople_CC))