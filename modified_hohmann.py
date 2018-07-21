#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np


def modified_hohmann(r1, r2):
    """
    Modified Hohmann Transfer
    """
    global mu, J2, req
    v1 = np.sqrt((mu/r1)*(1+(3*J2*req**2)/(2*r1**2)))
    v2 = np.sqrt((mu/r2)*(1+(3*J2*req**2)/(2*r2**2)))
    v_esc1 = np.sqrt((2*mu/r1)*(1+(J2*req**2)/(2*r1**2)))
    v_esc2 = np.sqrt((2*mu/r2)*(1+(J2*req**2)/(2*r2**2)))
    vA = np.sqrt(r2**2*(v_esc1**2-v_esc2**2)/(r2**2-r1**2))
    vB = np.sqrt(r1**2*(v_esc1**2-v_esc2**2)/(r2**2-r1**2))
    dV1 = vA - v1
    dV2 = v2 - vB
    return dV1, dV2, abs(dV1)+abs(dV2)


req = 6378
J2 = 1.08263e-3
mu = 398600
r1 = 6800
r2 = 8000
print(modified_hohmann(r1, r2))