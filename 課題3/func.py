import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
import math
import csv


def T_insert(T):
    rho_data = {460: 0.7667, 480: 0.7347,500: 0.7053, 550: 0.6412, 600: 0.5878}
    lambda_data = {460: 36.97, 480: 38.25, 500: 39.51, 550: 42.6, 600: 45.6}
    nu_data = {460: 33.51, 480: 36.01, 500: 38.58, 550: 45.27, 600: 45.6}
    Pr_data = {460: 0.711, 480: 0.710, 500: 0.710, 550: 0.709, 600: 0.710}
    if 460 <= T < 480:
        rho_insert = rho_data[460]*(480-T)/20+rho_data[480]*(T-460)/20
        lambda_insert = lambda_data[460]*(480-T)/20+lambda_data[480]*(T-460)/20
        nu_insert = nu_data[460]*(480-T)/20+nu_data[480]*(T-460)/20
        Pr_insert = Pr_data[460]*(480-T)/20+Pr_data[480]*(T-460)/20
        return rho_insert, lambda_insert, nu_insert, Pr_insert
    if 480 <= T < 500:
        rho_insert = rho_data[480]*(500-T)/20+rho_data[500]*(T-480)/20
        lambda_insert = lambda_data[480]*(500-T)/20+lambda_data[500]*(T-480)/20
        nu_insert = nu_data[480]*(500-T)/20+nu_data[500]*(T-480)/20
        Pr_insert = Pr_data[480]*(500-T)/20+Pr_data[500]*(T-460)/20
        return rho_insert, lambda_insert, nu_insert, Pr_insert
    if 500 <= T < 550:
        rho_insert = rho_data[500]*(550-T)/50+rho_data[550]*(T-500)/50
        lambda_insert = lambda_data[500]*(550-T)/50+lambda_data[550]*(T-500)/50
        nu_insert = nu_data[500]*(550-T)/50+nu_data[550]*(T-500)/50
        Pr_insert = Pr_data[500]*(550-T)/50+Pr_data[550]*(T-500)/50
        return rho_insert, lambda_insert, nu_insert, Pr_insert
    if 550 <= T < 600:
        rho_insert = rho_data[550]*(600-T)/50+rho_data[600]*(T-550)/50
        lambda_insert = lambda_data[550]*(600-T)/50+lambda_data[600]*(T-550)/50
        nu_insert = nu_data[550]*(600-T)/50+nu_data[600]*(T-550)/50
        Pr_insert = Pr_data[550]*(600-T)/50+Pr_data[600]*(T-550)/50
        return rho_insert, lambda_insert, nu_insert, Pr_insert



def heat_transfer(Q_flow,rho_T,d_haikan,Pr_T,L,nu_T,lambda_T,T_inf,T_w,S,T_C,C,eta_h,eta_m):
    """
    計算ここから
    """
    #T_m = (T_inf+T_w+273)/2
    #rho_T, lambda_T, nu_T, Pr_T = T_insert(T_m)
    u = Q_flow / rho_T / math.pi*4 / d_haikan / d_haikan
    Re = u*L / nu_T*(10**3)
    Nu_flatplate = 0.664*(Re**0.5)*(Pr_T**(1/3))
    h_transfer_flatplate = Nu_flatplate*lambda_T/L
    q_dot_flatplate = h_transfer_flatplate*(T_inf - T_w)
    Q_h_transfer_flatplate = q_dot_flatplate*S
    eta_net = (1-T_C/(T_w))*C*eta_h*eta_m
    L_net_flatplate = eta_net*Q_h_transfer_flatplate
    return Q_h_transfer_flatplate,eta_net,L_net_flatplate
