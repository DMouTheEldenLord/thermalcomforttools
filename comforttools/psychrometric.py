import numpy as np


def cal_pws(t_a):
    c1 = -5.6745359e3
    c2 = 6.3925247
    c3 = -9.677843e-3
    c4 = 6.2215701e-7
    c5 = 2.0747825e-9
    c6 = -9.484024e-13
    c7 = 4.1635019
    c8 = -5.8002206e3
    c9 = 1.3914993
    c10 = -4.8640239e-2
    c11 = 4.1764768e-5
    c12 = -1.4452093e-8
    c13 = 6.5459673
    pws = 0 * t_a
    i1 = (t_a < 0)
    i2 = (t_a >= 0)
    t = t_a[i1] + 273.15
    pws[i1] = np.exp(c1 / t + c2 + c3 * t + c4 * t ** 2 + c5 * t ** 3 + c6 * t ** 4 + c7 * np.log(t))
    t = t_a[i2] + 273.15
    pws[i2] = np.exp(c8 / t + c9 + c10 * t + c11 * t ** 2 + c12 * t ** 3 + c13 * np.log(t))
    return pws


def cal_h(t_a, d):
    h = 1.006 * t_a + 0.001 * d * (2501 + 1.86 * t_a)
    return h


def cal_d(p_a, b=101325):
    d = 621.945*p_a/(b - p_a)
    return d
