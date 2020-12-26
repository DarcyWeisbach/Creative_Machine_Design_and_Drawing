import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
import math
import csv

def task2(D_p,S_p,L_net,P_mean):
    V_SE = math.pi/4*D_p*D_p*S_p
    n = L_net/2.5/10**(-6)/P_mean/V_SE*1000
    return n,V_SE

def task3(delta_p,X_DE,V_SE,D_p):
    V_DE = X_DE*V_SE
    D_c = D_p+2*delta_p
    L_h = (V_DE-np.pi/4*D_c*D_c*2.5-1/3*np.pi/4*D_c*D_c*D_c /2*np.tan(np.pi/6))/(np.pi/4*(D_c*D_c-D_p*D_p))+2.5
    return L_h,D_c

def task4(d_k,X_DC,V_SE):
    V_DC = X_DC*V_SE
    L_k = V_DC/np.pi*4/d_k/d_k
    return L_k

def task5(S_p):
    L_s = 3*S_p
    return L_s

def task6(L_s,L_h):
    zeta = 1
    delta_h = 2.5
    delta_c = 1.5
    t_j = 8
    t_p = 0.5
    L_pe = zeta+L_s+t_p+t_j+t_p+L_h-delta_h
    L_pc = zeta+L_s+t_p-delta_c
    return L_pe,L_pc

def task7(D_c,D_p,V_SE):
    delta_c = 1.5
    t_j = 8
    V_R1 = np.pi/4*(D_c*D_c - D_p*D_p)*t_j
    V_R2 = np.pi/4*D_p*D_p*delta_c
    X_R = (V_R1+V_R2) / V_SE
    return X_R

def task8(X_DE,X_DC,X_R):
    X = X_DE+X_DC+X_R
    return X

def task9(T_w,T_C,kappa,alpha,P_mean,X):
    tau = T_C / T_w
    B = math.sqrt(tau*tau+2*tau*kappa*math.cos(alpha*math.pi/180)+kappa*kappa)
    S = tau+4*tau*X/(1+tau)+kappa
    delta = B/S
    P_max = P_mean*math.sqrt((1+delta)/(1-delta))
    return P_max,tau,delta

def task10(kappa,tau):
    phi = 180/math.pi*math.atan(kappa/tau)
    return phi

def task11(P_max,V_SE,delta,tau,n,eta_m,phi):
    W_i = P_max*(10**(-6))*V_SE*math.pi*delta*(1-tau)*math.sin(math.pi/180*phi) *math.sqrt((1 - delta)/(1 + delta))/(1+math.sqrt(1-delta*delta))
    L_i = W_i*n/60
    L_net_cal = L_i*eta_m
    return L_net_cal


