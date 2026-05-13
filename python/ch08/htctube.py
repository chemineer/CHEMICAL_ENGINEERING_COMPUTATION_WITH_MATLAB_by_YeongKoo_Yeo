def htctube(Nre, Pr, D, L, Xk, Phi):
    """
    Calculate convective heat transfer coefficient within tube (Python version of htctube.m)
    
    Args:
        Nre: Reynolds number
        Pr: Prandtl number
        D: inside diameter of tube (mm)
        L: tube length (m)
        Xk: thermal conductivity of fluid (W/m/K)
        Phi: viscosity correction factor
        
    Returns:
        Htube: convective heat transfer coefficient within tube
    """
    # mm -> m 단위 변환
    Dm = D * 1e-3 
    
    if Nre <= 2100:  # 층류(Laminar flow) 영역
        Gw = Nre * Pr * Dm / L
        if Gw > 100:
            Nu = 1.86 * Phi * (Gw**0.333)
        else:
            Nu = 3.66 + (0.085 * Gw * Phi) / (1 + 0.047 * (Gw**0.6667))
            
    elif Nre < 1e4:  # 전이(Transition) 영역
        Nu = 0.116 * (Nre**0.6667 - 125) * (Pr**0.333) * (1 + (Dm / L)**0.6667) * Phi
        
    else:  # 난류(Turbulent) 영역
        Nu = 0.023 * Phi * (Nre**0.8) * (Pr**0.333)
        
    # 열전달 계수 계산
    Htube = Nu * Xk / Dm
    
    return Htube