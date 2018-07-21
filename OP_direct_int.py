#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
numerically integration of Two Body Orbit Equation
"""
import numpy as np
from scipy.integrate import solve_ivp
import numpy.linalg as la
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
def dy_dt(t, f):
    global mu, J2, R
    x = f[0]
    y = f[1]
    z = f[2]
    vx = f[3]
    vy = f[4]
    vz = f[5]
    r = la.norm([x, y, z])
    px = (1.5*J2*mu*R**2/r**4)*(x/r)*(5*(z**2/r**2)-1)
    py = (1.5*J2*mu*R**2/r**4)*(y/r)*(5*(z**2/r**2)-1)
    pz = (1.5*J2*mu*R**2/r**4)*(z/r)*(5*(z**2/r**2)-3)
    ax = -mu*x/(r**3) + px
    ay = -mu*y/(r**3) + py
    az = -mu*z/(r**3) + pz
    dydt = [vx, vy, vz, ax, ay, az]
    return dydt

R = 6378
J2 = 0.00108263
t0 = 0
tf = 100*3600
mu = 398600
r0 = [1600, 5310, 3800]
v0 = [-7.35, 0.46, 2.47]
t = [t0, tf]
y0 = [r0[0], r0[1], r0[2], v0[0], v0[1], v0[2]]
tval = np.linspace(t0, tf, 10000)
sol = solve_ivp(dy_dt, t, y0, method='RK45', t_eval=tval,
                rtol=1e-12, atol=1e-12)
y = sol.y
y = y.tolist()
t = sol.t
t = t.tolist()
rx = y[0]
ry = y[1]
rz = y[2]
markers_on = rx[0]
descriptions = ["initial"]
ax = plt.axes(projection='3d')
ax.plot(rx, ry, rz, '-b.', markevery=markers_on)
plt.autoscale(enable=True, tight=True)
