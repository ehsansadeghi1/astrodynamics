from numpy import sin, cos, matrix
from constants import mu


def sv_from_coe(coe):
    """
    Calculation of the state vector from the orbital elements
    mu  - gravitational parameter (km^3/s^2)
    coe - orbital elements [h e RAAN i w TA]
        h    = angular momentum (km^2/s)
        e    = eccentricity
        RAAN = right ascension of the ascending node (rad)
        i = orbit inclination (rad)
        w    = argument of perigee (rad)
        TA   = true anomaly (rad)
    r   - position vector (km) in geocentric equatorial frame
    v   - velocity vector (km) in geocentric equatorial frame
    """
    h = coe[0]
    e = coe[1]
    RA = coe[2]
    i = coe[3]
    w = coe[4]
    TA = coe[5]
    dummy1 = matrix([[1], [0], [0]])
    dummy2 = matrix([[0], [1], [0]])
    rp = (h**2/mu)*(1/(1+e*cos(TA)))*(cos(TA)*dummy1+sin(TA)*dummy2)
    vp = (mu/h)*(-sin(TA)*dummy1+(e+cos(TA))*dummy2)
    del dummy1, dummy2

    R3_W = matrix([[cos(RA), sin(RA), 0],
                   [-sin(RA), cos(RA), 0],
                   [0, 0, 1]])

    R1_i = matrix([[1, 0, 0],
                   [0, cos(i), sin(i)],
                   [0, -sin(i), cos(i)]])

    R3_w = matrix([[cos(w), sin(w), 0],
                   [-sin(w), cos(w), 0],
                   [0, 0, 1]])
    Q_pX = (R3_w*R1_i*R3_W)
    Q_pX = Q_pX.transpose()
    r = Q_pX*rp
    r = [float(r[0]), float(r[1]), float(r[2])]
    v = Q_pX*vp
    v = [float(v[0]), float(v[1]), float(v[2])]
    return r, v
