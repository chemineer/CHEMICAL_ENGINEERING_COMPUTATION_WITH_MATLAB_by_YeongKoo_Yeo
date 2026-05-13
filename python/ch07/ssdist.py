import numpy as np

def ssdist(x, dpar):
    """
    증류탑 정상 상태 모델 (Steady-state Distillation Model)
    
    Args:
        x: 조성(mole fraction)을 담은 배열 (n개의 원소)
        dpar: 시스템 파라미터 딕셔너리
    Returns:
        f: 각 방정식의 잔차(residual) 배열
    """
    alpha = dpar['alpha']
    n = dpar['n']
    nf = dpar['nf']
    F = dpar['F']
    zf = dpar['zf']
    q = dpar['q']
    R = dpar['R']
    D = dpar['D']
    
    # 내부 유량 계산
    Lr = R
    B = F - D
    Ls = R + F * q
    Vs = Ls - B
    Vr = Vs + F * (1 - q)
    
    # 기상-액상 평형 관계 (y는 각 단의 기상 조성)
    y = alpha * x / (1 + (alpha - 1) * x)
    
    f = np.zeros(n)
    
    # 방정식 구성 (MATLAB의 1-based 인덱스를 파이썬 0-based로 변환)
    # 0번 인덱스: 응축기(Condenser)
    f[0] = Vr * y[1] - (D + R) * x[0]
    
    # 1번부터 nf-2번 인덱스: 정류부(Rectifying section)
    for i in range(1, nf - 1):
        f[i] = Lr * x[i-1] + Vr * y[i+1] - Lr * x[i] - Vr * y[i]
        
    # nf-1번 인덱스: 공급단(Feed tray)
    f[nf-1] = Lr * x[nf-2] + Vs * y[nf] - Ls * x[nf-1] - Vr * y[nf-1] + F * zf
    
    # nf부터 n-2번 인덱스: 회수부(Stripping section)
    for i in range(nf, n - 1):
        f[i] = Ls * x[i-1] + Vs * y[i+1] - Ls * x[i] - Vs * y[i]
        
    # n-1번 인덱스: 재비기(Reboiler)
    f[n-1] = Ls * x[n-2] - B * x[n-1] - Vs * y[n-1]
    
    return f