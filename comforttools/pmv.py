import math


def cal_pmv(ta=None, tr=None, rh=None, v=None, clo=None, met=None, top=None, w=0):
    """
        This function is used to calculated Predicted Mean Vote(PMV) through environmental parameters
        cal_pmv(ta=None, tr=None, rh=None, v=None, clo=None, met=None, top=None, w=0)
        ta: Air temperature, °C
        tr: Mean radiation temperature, °C
        rh: Relative humidity, %
        v:  Relative air velocity, m/s
        clo: Clothing insulation, clo
        met: Metabolic rate, met
        top: operative temperature, °C
        When top is NOT None, it will replace values of ta AND tr
        w: External work, W/m2

        example 1:
        Out[1]: cal_pmv(ta=25, tr=25, rh=50, v=0.1, clo=0.5, met=1.0)
        Out[1]: -0.4025360559202205

        example 2:
        Out[1]: cal_pmv(to=25, rh=50, v=0.1, clo=0.5, met=1.0)
        Out[1]: -0.4025360559202205
    """
    if top is not None:
        ta = top
        tr = top
    m = 58.15 * met
    if clo < 0.5:
        fcl = 1 + 0.2 * clo
    else:
        fcl = 1.05 + .1 * clo
    icl = 0.155 * clo
    tsk = 35.7 - 0.028 * m
    pa = rh / 100 * math.exp(16.6536 - 4030.183 / (ta + 235))
    tcl = min(ta, tsk)
    hcf = 12.1 * v ** 0.5
    max_iteration = int(1e5)
    hc = hcf
    hr = 4
    for i in range(0, max_iteration):
        hcn = 2.38 * abs(ta - tcl) ** 0.25
        hc = max(hcf, hcn)
        tcl_k = tcl + 273
        tr_k = tr + 273
        hr = 3.96e-8 * (tcl_k ** 2 + tr_k ** 2) * (tcl_k + tr_k)
        tcl1 = (icl * fcl * (hc * ta + hr * tr) + tsk) / (1 + icl * fcl * (hr + hc))
        if abs(tcl - tcl1) < 1e-4:
            break
        tcl = tcl1

    r = fcl * hr * (tcl - tr)
    c = fcl * hc * (tcl - ta)
    c_res = 0.0014 * m * (34 - ta)
    e_res = 0.017 * m * (5.867 - pa)
    e_dif = 3.05 * (5.733 - 0.00699 * (m - w) - pa)
    e_rsw = max(0, 0.42 * (m - w - 58.15))
    tl = m - w - c - r - e_rsw - e_res - e_dif - c_res
    pmv = (0.303 * math.exp(-0.036 * m) + 0.028) * tl
    return pmv


def cal_ppd(pmv):
    """
            This function is used to calculated Predicted Percentage of Dissatisfaction(PPD) through PMV
            PPD = 100 - 95 * exp(-(0.03353 * PMV^4 + 0.2179 * PMV^2))
    """
    return 100 - 95 * math.exp(-(0.03353 * pmv ** 4 + 0.2179 * pmv ** 2))


def cal_to_by_pmv(pmv=None, rh=None, v=None, clo=None, met=None, w=0):
    """
            This function is used to calculate operative temperature through PMV and other environmental parameters
            using Newton's method.

            example 1:
            Out[1]: cal_to(pmv=0, rh=50, v=0.1, clo=1.0, met=1.0)
            Out[1]: 23.28197017301155
    """
    max_iteration = int(1e5)
    top = 25
    dt = 0.1
    flag = 0
    for i in range(0, max_iteration):
        f1 = cal_pmv(rh=rh, v=v, clo=clo, met=met, top=top, w=w) - pmv
        f2 = cal_pmv(rh=rh, v=v, clo=clo, met=met, top=top+dt, w=w) - pmv
        dx = f1/(f2-f1) * dt
        if abs(dx) < 1e-3:
            flag = 1
            break
        top = top - dx
    if flag == 0:
        return -999
    else:
        return top
