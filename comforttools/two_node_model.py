import math

met_factor = 58.2
sbc = 5.6697e-8


class TwoNodeModel:
    h_sk = 0
    w = 0
    ps_sk = 0
    hr = 4.7
    met = 0
    work = 0

    def __init__(self, height=1.8, weight=70,
                 c_sw=170, c_dil=120, c_str=0.5,
                 t_sk_0=33.7, t_cr_0=36.8, t_b_0=36.49,
                 skin_blood_flow_0=6.3, alpha_0=0.1):
        self.height = height
        self.weight = weight
        self.a_du = 0.20247 * self.height ** 0.725 * self.weight ** 0.425
        self.c_sw = c_sw
        self.c_dil = c_dil
        self.c_str = c_str
        self.t_sk_0 = t_sk_0
        self.t_sk = t_sk_0
        self.t_cr_0 = t_cr_0
        self.t_cr = t_cr_0
        self.t_b_0 = t_b_0
        self.t_b = t_b_0
        self.skin_blood_flow_0 = skin_blood_flow_0
        self.skin_blood_flow = skin_blood_flow_0
        self.alpha_0 = alpha_0
        self.alpha = alpha_0

    def two_node_set(self, t_sk, t_cr,
                     skin_blood_flow, alpha):
        self.t_sk = t_sk
        self.t_cr = t_cr
        self.skin_blood_flow = skin_blood_flow
        self.alpha = alpha

    def reset(self):
        self.two_node_set(t_sk=self.t_sk_0, t_cr=self.t_cr_0,
                          skin_blood_flow=self.skin_blood_flow_0,
                          alpha=self.alpha_0)

    def simulate(self, ta=None, tr=None, v=None, rh=None, met=None, clo=None, b=101.325, top=None, icl=0.45,
                 work=0, duration_time=100, ashrae=True):
        if ashrae is True:
            self.weight = 69.9
            self.a_du = 1.8258
            duration_time = 60
        if top is not None:
            ta = top
            tr = top
        b = b*0.009869
        pw = rh/100*cal_pws(ta)
        v = max(v, 0.1)
        lr = 2.2/b
        rcl = clo*0.155
        fcl = 1 + 0.15*clo
        self.met = met
        self.work = work
        rm = met*met_factor
        m = met*met_factor
        if clo <= 0:
            w_crit = 0.38*v**(-0.29)
            icl = 1
        else:
            w_crit = 0.59*v**(-0.08)
        hc_n = 3*b**0.53
        hc_c = 8.600001*(v*b)**0.53
        hc = max(hc_n, hc_c)
        ht = hc + self.hr
        ra = 1/(fcl*ht)   # 空气层热阻
        top = (self.hr*tr + hc*ta)/ht
        tcl = top + (self.t_sk - top)/(ht*(ra+rcl))  # (tcl-top)*ht = (t_sk - top)/(ra+rcl) 似乎不太对
        e_sk = met * 0.1
        dry = 0
        p_wet = 0
        flag = 1
        tcl1 = tcl
        for i in range(0, duration_time):
            # 传热计算
            while 1:    # 迭代确定辐射换热系数
                if flag == 1:
                    tcl1 = tcl
                    self.hr = 4 * sbc * ((tcl + tr) / 2 + 273.15) ** 3 * 0.72
                    ht = hc + self.hr
                    ra = 1 / (fcl * ht)
                    top = (self.hr * tr + hc * ta) / ht
                tcl = (ra*self.t_sk + rcl*top)/(ra + rcl)    # (tcl-top)/ra = (t_sk-tcl)/rcl
                flag = 1
                if abs(tcl-tcl1) <= 1e-2:
                    break
            flag = 0
            dry = (self.t_sk - top)/(ra+rcl)
            hfcs = (self.t_cr - self.t_sk)*(5.28 + 1.163*self.skin_blood_flow)     # 血流换热量
            e_res = 0.0023*m*(44-pw)
            c_res = 0.0014*m*(34-ta)
            s_cr = m - hfcs - e_res - c_res - work
            s_sk = hfcs - dry - e_sk
            tc_sk = 0.97*self.alpha*self.weight
            tc_cr = 0.97*(1-self.alpha)*self.weight
            dt_sk = (s_sk*self.a_du)/(tc_sk*60)
            dt_cr = (s_cr*self.a_du)/(tc_cr*60)
            self.t_sk = self.t_sk + dt_sk
            self.t_cr = self.t_cr + dt_cr
            self.t_b = self.alpha*self.t_sk + (1 - self.alpha)*self.t_cr
            # 控制信号计算
            sk_sig = self.t_sk - self.t_sk_0
            warm_s = max(0.0, sk_sig)
            cold_s = max(0.0, -sk_sig)
            cr_sig = self.t_cr - self.t_cr_0
            warm_c = max(0.0, cr_sig)
            cold_c = max(0.0, -cr_sig)
            bd_sig = self.t_b-self.t_b_0
            warm_b = max(0.0, bd_sig)
            self.skin_blood_flow = (self.skin_blood_flow_0 + self.c_dil*warm_c)/(1 + self.c_str*cold_s)
            self.skin_blood_flow = max(0.5, min(90.0, self.skin_blood_flow))
            reg_sw = self.c_sw*warm_b*math.exp(warm_s/10.7)
            reg_sw = min(reg_sw, 500.0)
            e_rsw = 0.68*reg_sw
            rea = 1/(lr*fcl*hc)
            recl = rcl/(lr*icl)
            e_max = (cal_pws(self.t_sk) - pw)/(rea + recl)
            p_rsw = e_rsw/e_max
            p_wet = 0.06 + 0.94*p_rsw
            e_dif = p_wet*e_max - e_rsw
            if p_wet > w_crit:
                p_wet = w_crit
                p_rsw = w_crit/0.94
                e_rsw = p_rsw*e_max
                e_dif = 0.06*(1-p_rsw)*e_max
            if e_max < 0:
                e_dif = 0
                e_rsw = 0
                p_wet = w_crit
            e_sk = e_rsw + e_dif
            m_shiv = 19.4*cold_s*cold_c
            m = rm + m_shiv
            self.alpha = 0.0417737+0.7451833/(self.skin_blood_flow+0.585417)
        self.h_sk = dry + e_sk
        self.w = p_wet
        self.ps_sk = cal_pws(self.t_sk)


def cal_pws(t):
    return math.exp(18.6686-4030.183/(t+235))


def cal_set(ta=None, tr=None, v=None, rh=None, met=None, clo=None, top=None, b=101.325, work=0):
    model = TwoNodeModel()
    model.simulate(ta=ta, tr=tr, v=v, rh=rh, met=met, clo=clo, top=top, b=b, work=work)
    b = b * 0.009869
    lr = 2.2 / b
    hr_s = model.hr
    if model.met < 0.85:
        hc_s = 3
    else:
        hc_s = max(5.66 * (model.met - 0.85) ** 0.39, 3)
    clo_s = 1.52 / ((model.met - work / met_factor) + 0.6944) - 0.1835
    ht_s = hc_s + hr_s
    rcl_s = 0.155 * clo_s
    facl_s = 1 + 0.25 * clo_s  #
    im_s = 0.45
    ra_s = 1 / (facl_s * ht_s)
    rea_s = 1 / (lr * facl_s * hc_s)
    fcl_s = 1 / (1 + 0.155 * facl_s * ht_s * clo_s)
    icl_s = im_s * hc_s / ht_s * (1 - fcl_s) / (hc_s / ht_s - fcl_s * im_s)
    recl_s = rcl_s / (lr * icl_s)
    hd_s = 1 / (ra_s + rcl_s)
    he_s = 1 / (rea_s + recl_s)
    # 用牛顿法求SET
    dt = 1e-4
    dx = 100
    set0 = model.t_sk - model.h_sk / hd_s
    while abs(dx) > 1e-2:
        err1 = model.h_sk - hd_s * (model.t_sk - set0) - model.w * he_s * (model.ps_sk - 0.5 * cal_pws(set0))
        err2 = model.h_sk - hd_s * (model.t_sk - (set0 + dt)) - model.w * he_s * (model.ps_sk - 0.5 * cal_pws(set0 + dt))
        dx = err1 / (err2 - err1) * dt
        set0 = set0 - dx
    return set0
