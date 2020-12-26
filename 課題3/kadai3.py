import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
import math
import func
import csv
import calc

"""
設計条件 ここから
"""
kappa = 1
X_DE = 1.5
X_DC = 0.5
P_mean = 129
#T_E = 673
#今回の条件では，T_wである．

T_C = 273+50
alpha = 90
#L_net = 1.2
#S_p = 8
#D_p = 10
D_p_list = [10, 12, 15, 18, 22, 25]
L_list = {10: 60, 12: 68, 15: 89, 18: 105, 22: 111, 25: 142}
C = 0.10
eta_h = 0.50
eta_m = 0.80
d = 2
d_haikan = 108.3*10**(-3)
L = 50
t = 50
S = (2*L*d+d*t+2*L*t)*(10**(-6))*8
Q_flow = 0.04
T_inf = 300+273
data=[["温度 [K]","L","t","d","面積S","計算したQ_h","L_net(設計条件となる)","D_p","S_p","delta_p","d_k","回転数","L_h","L_k","L_s","L_pe","L_pc","X_R","X","P_max","phi","L_net"]]

"""
設計条件ここまで
"""

"""
設計計算ここから
"""
for T_w in range(200,250,1):
    #温度ごとにする
    flg=True
    T_w=T_w+273
    T_m = (T_inf+T_w)/2
    rho_T, lambda_T, nu_T, Pr_T = func.T_insert(T_m)
    Q_h_transfer_flatplate, eta_net, L_net_flatplate = func.heat_transfer(Q_flow, rho_T, d_haikan, Pr_T, L, nu_T, lambda_T, T_inf, T_w, S, T_C, C, eta_h, eta_m)
    if L_net_flatplate<=1.2:
        continue
    for D_p in D_p_list:
        for S_p in np.arange(D_p-0.1,6,-0.1):
            S_p=round(S_p,1)
            if (S_p>=22.0):
                continue
            #n, V_SE = calc.task2(D_p, S_p, L_net, P_mean)
            n, V_SE = calc.task2(D_p, S_p, L_net_flatplate, P_mean)
            if (n>2000):
                continue
            for delta_p in np.arange(3.0,0.75,-0.1):
                delta_p=round(delta_p,2)
                L_h, D_c = calc.task3(delta_p, X_DE, V_SE, D_p)
                if (L_h > 80):
                    continue
                for d_k in np.arange(2.5,5.0,0.1):
                    d_k = round(d_k,1)
                    L_k = calc.task4(d_k, X_DC, V_SE)
                    if (L_k>100):
                        continue
                    L_s = calc.task5(S_p)
                    L_pe, L_pc = calc.task6(L_s, L_h)
                    if (L_pe > L_list[D_p] or L_pc > L_list[D_p]):
                        continue
                    X_R = calc.task7(D_c, D_p, V_SE)
                    X = calc.task8(X_DE, X_DC, X_R)
                    P_max, tau, delta = calc.task9(T_w, T_C, kappa, alpha, P_mean,X)
                    phi = calc.task10(kappa, tau)
                    L_net_cal = calc.task11(P_max, V_SE, delta, tau, n, eta_m,phi)
                    if (L_net_cal < 1.2):
                        continue
                    elif (L_net_cal >= 1.2):
                        if (flg==True):
                            data.append([T_w, L, t, d, S, Q_h_transfer_flatplate, L_net_flatplate, D_p, S_p,delta_p, d_k, n, L_h, L_k, L_s, L_pe, L_pc, X_R, X, P_max, phi, L_net_cal])
                            print("exist")
                            flg=False
                      
with open('8枚_温度ごと.csv', 'w') as file:
    writer = csv.writer(file, lineterminator='\n')
    writer.writerows(data)
print("finish")
