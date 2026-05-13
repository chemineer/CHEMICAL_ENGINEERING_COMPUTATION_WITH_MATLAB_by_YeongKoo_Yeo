import numpy as np
from types import SimpleNamespace

def rfobj(rf, opdat):
    """
    rfobj.m: 환류비(rf)에 따른 증류탑 비용 및 목적 함수 계산
    
    Inputs:
    rf    : 환류비 (Reflux ratio)
    opdat : 파라미터가 담긴 객체 (SimpleNamespace 등)
    
    Outputs:
    C     : 총 비용 (목적 함수 값)
    """
    # 데이터 추출
    we = opdat.we
    eta = opdat.eta
    S = opdat.S
    K = opdat.K
    rhoG = opdat.rhoG
    rhoL = opdat.rhoL
    rhos = opdat.rhos
    lamb = opdat.lamb
    lambs = opdat.lambs
    Css = opdat.Css
    Cst = opdat.Cst
    F = opdat.F
    T = opdat.T
    P = opdat.P
    alpa = opdat.alpa
    xB = opdat.xB
    xD = opdat.xD
    xF = opdat.xF
    
    # 파라미터 설정
    rLV = rf / (rf + 1)
    rDV = 1 / (1 + rf)
    rpLV = (rf + ((xD - xB) / (xF - xB))) / (1 + rf)
    rBV = (((xD - xB) / (xF - xB)) - 1) / (1 + rf)
    
    # 트레이 수 계산
    N = 1
    ye = [xD]
    xe = [ye[0] / (alpa * (1 - ye[0]) + ye[0])]
    
    while xe[N-1] > xB:
        if xe[N-1] > xF:
            # Rectifying section line
            new_ye = rLV * xe[N-1] + rDV * xD
        else:
            # Stripping section line
            new_ye = rpLV * xe[N-1] - rBV * xB
        
        ye.append(new_ye)
        xe.append(new_ye / (alpa * (1 - new_ye) + new_ye))
        N += 1
        
    # 탑 설계 계산
    V = F * (1 + rf) * (xF - xB) / (xD - xB) # 증기 속도
    d = np.sqrt((4 * V * 22.4 * 760 * (T + 273.15)) / (273 * np.pi * P * K * np.sqrt((rhoL - rhoG) / rhoL)))
    h = 0.6 * ((N - 1) / eta + 1) + 2       # 높이
    w = 14.7 * (P / 760) * d / 2 / (we * S - 0.6 * 14.7 * (P / 760)) + 0.0032 # 강철 두께
    A = 4 * np.pi * (d / 2)**2 + np.pi * d * h # 표면적
    
    # 비용 및 목적 함수 계산
    Ceng = lamb * F * (1 + rf) * Css * (xF - xB) / (xD - xB) / lambs # 에너지 비용
    Ccol = (A + np.pi * N * (d / 2)**2) * w * rhos * Cst            # 컬럼 비용
    C = Ccol / 3 + Ceng
    
    return C

# --- 사용 예시 ---
# opdat = SimpleNamespace(we=..., eta=..., S=..., K=..., ...)
# result = rfobj(2.5, opdat)