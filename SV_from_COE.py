#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from numpy import sin, cos, matrix, transpose, deg2rad


def sv_from_coe(coe, mu):
    """
    Calculation of the state vector from the orbital elements
    mu  - gravitational parameter (km^3/s^2)
    coe - orbital elements [a e RAAN i w TA]
    where
        h    = angular momentum (km^2/s)
        e    = eccentricity
        RAAN = right ascension of the ascending node (rad)
        incl = orbit inclination (rad)
        w    = argument of perigee (rad)
        TA   = true anomaly (rad)
    r   - position vector (km) in geocentric equatorial frame
    v   - velocity vector (km) in geocentric equatorial frame
    """
    h = coe[0]
    e = coe[1]
    RA = coe[2]
    incl = coe[3]
    w = coe[4]
    TA   = coe[5]
    dummy1 = matrix([[1], [0], [0]])
    dummy2 = matrix([[0], [1], [0]])
    rp = (h**2/mu)*(1/(1+e*cos(TA)))*(cos(TA)*dummy1+sin(TA)*dummy2);
    vp = (mu/h)*(-sin(TA)*dummy1+(e+cos(TA))*dummy2)

    R3_W = matrix([[cos(RA), sin(RA), 0],
                      [-sin(RA), cos(RA), 0],
                      [0, 0, 1]])

    R1_i = matrix([[1, 0, 0],
                      [0, cos(incl), sin(incl)],
                      [0, -sin(incl), cos(incl)]])

    R3_w = matrix([[cos(w), sin(w), 0],
                      [-sin(w), cos(w), 0],
                      [0, 0, 1]])
    Q_pX = (R3_w*R1_i*R3_W)
    Q_pX = Q_pX.transpose()
    r = Q_pX*rp
    v = Q_pX*vp
    r = [float(r[0]), float(r[1]), float(r[2])]
    v = [float(v[0]), float(v[1]), float(v[2])]
    return r, v
mu   = 398600
h    = 80000;
e    = 1.4;
RA   = 40;
incl = 30;
w    = 60;
TA   = 30;
coe = [h, e, deg2rad(RA), deg2rad(incl), deg2rad(w), deg2rad(TA)]
r = sv_from_coe(coe, mu)[0]
v = sv_from_coe(coe, mu)[1]
print(r)
print(v)