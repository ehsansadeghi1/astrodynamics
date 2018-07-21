#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np


def bi_elliptic(r_A, r_B, r_C):
    """
    bi-elliptic transfers
    """
    global mu
    a = r_C/r_A
    b = r_B/r_A
    v0 = np.sqrt(mu/r_A)
    dV_BE = np.sqrt(2*(a+b)/(a*b))-(1+np.sqrt(a))/np.sqrt(a)-np.sqrt(2/(b*(1+b)))*(1-b)
    dV = dV_BE*v0
    return dV

mu = 398600
r_A = 7000
r_C = 105000
r_B = 210000
print(bi_elliptic(r_A, r_B, r_C))