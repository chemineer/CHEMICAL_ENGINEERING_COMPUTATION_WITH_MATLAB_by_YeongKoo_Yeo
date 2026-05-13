import numpy as np
from scipy.optimize import fsolve
from satsteam import satsteam

def multievapSI(evdat):
    # Data
    crit = 1e-3
    Tc = 647.096
    Pc = 22064000
    T0 = 273.15
    xf = evdat['xf']
    xpn = evdat['xp']
    mf = evdat['mf']
    Tf = evdat['Tf'] + T0
    Ps = evdat['Ps']
    Pn = evdat['Pn']
    U = np.array(evdat['U'])
    n = len(U)
    
    # Saturation temperature calculation
    a = np.array([-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502])
    
    def eq(x, P):
        return (a[0]*x + a[1]*x**1.5 + a[2]*x**3 + a[3]*x**3.5 + 
                a[4]*x**4 + a[5]*x**7.5 - (1-x)*np.log(P/Pc))
    
    # fzero 대신 scipy.optimize.fsolve 사용
    xs = fsolve(eq, 0.5, args=(Ps,))[0]
    xn = fsolve(eq, 0.5, args=(Pn,))[0]
    
    Ts = Tc * (1 - xs)
    Tn = Tc * (1 - xn)
    
    # Enthalpy
    sts = satsteam(Ts - T0)
    stv = satsteam(Tn - T0)
    stf = satsteam(Tf - T0)
    
    Hv0, hp0 = sts['hV'], sts['hL']
    Hvn, hpn, hf = stv['hV'], stv['hL'], stf['hL']
    
    # Mass balances and temperature drop
    mpn = xf * mf / xpn
    dTtotal = Ts - Tn
    sumU = np.sum(1.0 / U)
    dT = (1.0 / U) * dTtotal / sumU
    
    critA = 10
    oldA = 10 * np.ones(n)
    iter_count = 0
    
    T = np.zeros(n)
    
    while critA >= crit:
        T[0] = Ts - dT[0]
        for j in range(1, n):
            T[j] = T[j-1] - dT[j]
            
        Hv = np.zeros(n)
        hp = np.zeros(n)
        for j in range(n):
            sprop = satsteam(T[j] - T0)
            Hv[j] = sprop['hV']
            hp[j] = sprop['hL']
            
        # balance equations (matrix system: evM * mv = evb)
        evM = np.array([
            [Hv0 - hp0, hp[0] - Hv[0], 0],
            [0, Hv[0] + hp[1] - 2*hp[0], hp[1] - Hv[1]],
            [0, Hv[2] - hp[1], Hv[1] + Hv[2] - 2*hp[1]]
        ])
        evb = np.array([hp[0] - hf, hp[1] - hp[0], Hv[2] - hp[1] + (hp[2] - Hv[2])*xf/xpn]) * mf
        
        mv = np.linalg.solve(evM, evb)
        
        q = np.zeros(n)
        q[0] = mv[0] * (Hv0 - hp0)
        for j in range(1, n):
            q[j] = mv[j] * (Hv[j] - hp[j])
            
        Area = q / (U * dT)
        avgA = np.sum(Area) / n
        critA = np.sum(np.abs(Area - oldA))
        dT = dT * Area / avgA
        
        if np.abs(np.sum(dT) - dTtotal) >= crit:
            dT = dT * dTtotal / np.sum(dT)
            
        iter_count += 1
        oldA = Area
        
    res = {
        'T': T - T0,
        'A': Area,
        'mv': mv,
        'iter': iter_count
    }
    return res