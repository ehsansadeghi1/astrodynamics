#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np


def Tsiolkovsky(I_sp, m_0, m_f):
    """
    Tsiolkovsky rocket equation
    m_0 = initial mass
    m_f = final mass
    """
    global g_0
    U_e = I_sp*g_0
    deltaV = U_e*np.log(m_0/m_f)
    return deltaV


def Hohmann_transfer(e1, r_p1, e2, r_p2):
    global mu
    r_a1 = r_p1*((1+e1)/(1-e1))
    r_a2 = r_p2*((1+e2)/(1-e2))
    """
    h = angular mpmentum
    1, 2: initial and final orbits
    3, 4: two different transfer orbits
    h = np.sqrt(2*mu)*np.sqrt(r_a*r_p/(r_a+r_p))
    """
    h1 = np.sqrt(2*mu)*np.sqrt(r_a1*r_p1/(r_a1+r_p1))
    h2 = np.sqrt(2*mu)*np.sqrt(r_a2*r_p2/(r_a2+r_p2))
    h3 = np.sqrt(2*mu)*np.sqrt(r_a2*r_p1/(r_a2+r_p1))
    h4 = np.sqrt(2*mu)*np.sqrt(r_a1*r_p2/(r_a1+r_p2))
    V_p1 = h1/r_p1
    V_a1 = h1/r_a1
    V_p2 = h2/r_a2
    V_a2 = h2/r_a2
    V_p3 = h3/r_p1
    V_a3 = h3/r_a2
    V_p4 = h4/r_a1
    V_a4 = h4/r_p2
    dV1 = abs(V_p1-V_p3)+abs(V_a2-V_a3)  # orbit 3
    dV2 = abs(V_a1-V_p4)+abs(V_p2-V_a4)  # orbit 4
    return min(dV1, dV2)
