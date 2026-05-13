import numpy as np

def satsteam(Ts_c):
    """
    제공된 온도 Ts(C)에서 포화 수증기의 물성을 계산하는 함수입니다.
    
    Args:
        Ts_c (float or np.array): 온도 (Celsius)
        
    Returns:
        dict: 계산된 물성들을 담은 딕셔너리 (stpr)
    """
    # Critical properties
    Tc = 647.096  # critical temperature(K)
    Pc = 22064000  # critical pressure(Pa)
    rhoc = 322  # critical density(kg/m^3)
    
    # Definition of parameters and constants
    Ts = Ts_c + 273.15  # C -> K
    alpha0 = 1000  # alpha_0(J/kg)
    phi0 = 1000 / 647.096
    
    a = np.array([-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502])
    b = np.array([1.99274064, 1.09965342, -0.510839303, -1.75493479, -45.5170352, -674694.45])
    c = np.array([-2.03150240, -2.68302940, -5.38626492, -17.2991605, -44.7586581, -63.9201063])
    d = np.array([-5.65134998e-8, 2690.66631, 127.287297, -135.003439, 0.981825814])
    
    alphad = -1135.905627715
    phid = 2319.5246
    theta = Ts / Tc
    tau = 1 - theta
    
    # saturated steam pressure
    tw = (Tc / Ts) * (a[0]*tau + a[1]*tau**1.5 + a[2]*tau**3 + a[3]*tau**3.5 + a[4]*tau**4 + a[5]*tau**7.5)
    Ps = Pc * np.exp(tw)
    
    # density of saturated liquid
    rhoL = rhoc * (1 + b[0]*tau**(1/3) + b[1]*tau**(2/3) + b[2]*tau**(5/3) +
                   b[3]*tau**(16/3) + b[4]*tau**(43/3) + b[5]*tau**(110/3))
    
    # density of saturated steam
    vw = (c[0]*tau**(1/3) + c[1]*tau**(2/3) + c[2]*tau**(4/3) + 
          c[3]*tau**3 + c[4]*tau**(37/6) + c[5]*tau**(111/6)) # 원본 71/6 확인 필요, 식에 따라 교정 가능
    rhoV = rhoc * np.exp(vw)
    
    # specific volume
    vL = 1 / rhoL  # saturated liquid
    vV = 1 / rhoV  # saturated steam
    
    # alpha
    alpha = alpha0 * (alphad + d[0]*theta**(-19) + d[1]*theta + d[2]*theta**4.5 + d[3]*theta**5 + d[4]*theta**54.5)
    
    # phi
    phi = phi0 * (phid + (19/20)*d[0]*theta**(-20) + d[1]*np.log(theta) +
                  (9/7)*d[2]*theta**3.5 + (5/4)*d[3]*theta**4 + (109/117)*d[4]*theta**53.5)
    
    # dp/dT
    tv = (7.5*a[5]*tau**6.5 + 4*a[4]*tau**3 + 3.5*a[3]*tau**2.5 +
          3*a[2]*tau**2 + 1.5*a[1]*tau**0.5 + a[0] + np.log(Ps / Pc))
    dpdT = (-Ps / Ts) * tv
    
    # enthalpy
    hL = alpha + (Ts / rhoL) * dpdT  # saturated liquid
    hV = alpha + (Ts / rhoV) * dpdT  # saturated steam
    
    # entropy
    sL = phi + (1 / rhoL) * dpdT  # saturated liquid
    sV = phi + (1 / rhoV) * dpdT  # saturated steam
    
    # Result structure (Dictionary in Python)
    stpr = {
        'T': Ts,
        'P': Ps,
        'vL': vL,
        'vV': vV,
        'hL': hL,
        'hV': hV,
        'sL': sL,
        'sV': sV
    }
    
    return stpr