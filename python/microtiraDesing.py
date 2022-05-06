
import numpy as np
import cmath
def microtira(tipo,t,H,er,zo):
    if (tipo):
        # Hammerstad
        A = (zo/60 * np.sqrt((er+1)/2) ) + (((er-1)/(er+1)) * (0.23 + (0.11/er)))
        B = ((377 * np.pi) / (2*zo*np.sqrt(er)))    
        print("Hammerstad".capitalize().center(50,"-"))
        print('A=',A)
        print('B=',B)

    else:
        # Wheeler
        A =  (zo/60 * np.sqrt((er+1)/2) ) + (((er-1)/(er+1)) * (0.226 + (0.121/er)))
        B = ((377 * np.pi) / (2*zo*np.sqrt(er)))
        print("Wheeler".capitalize().center(50,"-"))
        print('A=',A)
        print('B=',B)


    if(tipo):
        # Hammerstad
        W_H_ratio_1 = (8*np.exp(A))/(np.exp(2*A)-2)
        W_H_ratio_2 = (2/np.pi) * (B-1 - np.log((2*B)-1)+((er-1)/(2*er)) * (np.log(B-1)+0.39-0.61/er))

    else:
        # Wheeler
        W_H_ratio_1 = (8*np.exp(A))/(np.exp(2*A)-2)
        W_H_ratio_2 = (((er-1)/(np.pi*er))* (np.log(B-1) + 0.293 - (0.517/er))+ (2/np.pi) *(B-1-np.log(2*B-1)))

    if (W_H_ratio_1 <= 2):
        # print("Relacion W/H menor o igual a 2")
        W = H * W_H_ratio_1
    if (W_H_ratio_2 >= 2):
        # print("Relacion W/H mayor o igual a 2")
        W = H * W_H_ratio_2

    print('W[mm]= ',W)

    ratio = W/H
    print('W/H = ',ratio)
    if(ratio <= (1/(2*np.pi))):
        # print('Relacion W/H <1/2*pi')
        We = W + ((t/np.pi)*(1+np.log((4*np.pi*W)/t)))
    if(ratio >= (1/(2*np.pi))):
        # print('Relacion W/H >1/2*pi')
        We = W + ((t/np.pi)*(1+np.log((2*H)/t)))
    print('We[mm]= ',We)

    if(ratio<=1):
        # print('Relacion W/H <1')
        er_p = (er+1)/2 + (((er-1)/2)*((1/(np.sqrt(1+(12*H/W))))+0.04*(1-W/H)**2))
    if(ratio>=1):
        # print('Relacion W/H >1')
        er_p = (er+1)/2 + (((er-1)/2)*(1/(np.sqrt(1+(12*H/W)))))
    print('er_p= ',er_p)

    if(ratio<=1):
        zo_p = 60/np.sqrt(er_p) * np.log(8*H/W + W/4*H)
    if(ratio>=1):
        if(tipo):
            zo_p = (120*np.pi/np.sqrt(er_p)) /  (W/H + 1.393 + 0.667*np.log(W/H + 1.444))
        else: 
            zo_p = (120*np.pi/np.sqrt(er_p)) /  (W/H + 2.46 - 0.49*H/W + (1-H/W)**6)

    print('Impedancia obtenida: ',zo_p)

    return A,B,W,ratio,er_p,zo_p
