import numpy as np

def dyndist(t, x, dpar, dels):
    """
    이성분 증류탑을 위한 미분 방정식 시스템 구현
    t: 현재 시간
    x: 액체 몰 분율 벡터 (상태 변수)
    dpar: 기본 운전 매개변수 (딕셔너리)
    dels: 운전 조건 변화 매개변수 (딕셔너리)
    """
    # 매개변수 추출
    alpha = dpar['alpha']
    n = dpar['n']
    nf = dpar['nf']
    Fi = dpar['F']
    zfi = dpar['zf']
    q = dpar['q']
    Ri = dpar['R']
    Vsi = dpar['Vs']
    md = dpar['md']
    mb = dpar['mb']
    mt = dpar['mt']
    
    delR = dels['delR']; delRt = dels['delRt']
    delV = dels['delV']; delVt = dels['delVt']
    delz = dels['delz']; delzt = dels['delzt']
    delF = dels['delF']; delFt = dels['delFt']
    
    # 시간(t)에 따른 운전 조건 변화 적용 (Step change)
    R = Ri + delR if t >= delRt else Ri
    Vs = Vsi + delV if t >= delVt else Vsi
    zf = zfi + delz if t >= delzt else zfi
    F = Fi + delF if t >= delFt else Fi
    
    # 유량 계산
    Lr = R
    Ls = R + F * q
    B = Ls - Vs
    D = F - B
    Vr = Vs + F * (1 - q)
    
    # 초기화 및 상평형 관계 (y = alpha*x / (1 + (alpha-1)*x))
    dx = np.zeros(n)
    y = alpha * x / (1 + (alpha - 1) * x)
    
    # 파이썬 인덱스는 0부터 시작하므로 MATLAB 인덱스 i를 i-1로 조정하여 적용
    
    # 1. 응축기 (Condenser, MATLAB x(1))
    dx[0] = Vr * (y[1] - x[0]) / md
    
    # 2. 농축부 (Rectifying section, MATLAB 2:nf-1)
    # nf는 MATLAB 기준이므로 파이썬 범위는 range(1, nf-1)
    for i in range(1, nf - 1):
        dx[i] = (Lr * x[i-1] + Vr * y[i+1] - Lr * x[i] - Vr * y[i]) / mt
        
    # 3. 원료 주입단 (Feed stage, MATLAB x(nf))
    # MATLAB 인덱스 nf는 파이썬 인덱스 nf-1에 대응
    dx[nf-1] = (Lr * x[nf-2] + Vs * y[nf] + F * zf - Ls * x[nf-1] - Vr * y[nf-1]) / mt
    
    # 4. 회수부 (Stripping section, MATLAB nf+1:n-1)
    # 파이썬 범위는 range(nf, n-1)
    for i in range(nf, n - 1):
        dx[i] = (Ls * x[i-1] + Vs * y[i+1] - Ls * x[i] - Vs * y[i]) / mt
        
    # 5. 재비기 (Reboiler, MATLAB x(n))
    # MATLAB 인덱스 n은 파이썬 인덱스 n-1에 대응
    dx[n-1] = (Ls * x[n-2] - B * x[n-1] - Vs * y[n-1]) / mb
    
    return dx