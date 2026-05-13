import numpy as np

def gaopt(fun, xl, xu, nb, ps, ng, mp):
    """
    gaopt: maximization by the genetic algorithm
    """
    n = len(xu)
    # 초기 인구 생성 (0과 1의 랜덤 행렬)
    Pn = np.random.randint(2, size=(ps, n * nb))
    
    xbest = np.zeros(n)
    fbest = -np.inf
    
    for kg in range(ng):
        ft = fitness(fun, ps, nb, n, xl, xu, Pn)
        fmax_idx = np.argmax(ft)
        fmax = ft[fmax_idx]
        
        # 최적값 업데이트
        if kg == 0 or fmax > fbest:
            fbest = fmax
            # 최적해 디코딩
            for k in range(n):
                bits = Pn[fmax_idx, nb*k : nb*(k+1)]
                adec = int("".join(map(str, bits)), 2)
                xbest[k] = xl[k] + (xu[k] - xl[k]) / (2**nb - 1) * adec
        
        # 스케일링 및 룰렛 휠 계산
        frk, ft_sorted = frank(ps, ft)
        am = np.zeros(ps)
        for k in range(ps):
            ik = int(frk[k]) - 1 # 0-indexed
            am[ik] = 2 * (ps + 1 - (k + 1)) / (ps * (ps + 1))
        
        # 누적 확률 계산
        am = np.cumsum(am)
        
        # 선택 및 교차/변이
        new_Pn = np.zeros_like(Pn)
        for k in range(ps):
            kstr1 = roul(ps, am)
            kstr2 = roul(ps, am)
            kchd = cros(Pn, nb, n, kstr1, kstr2)
            kchd = mutat(kchd, mp, nb, n)
            new_Pn[k, :] = kchd
        Pn = new_Pn
        
    return xbest, fbest, ng * ps

def fitness(fun, ps, nb, n, xl, xu, Pn):
    ft = np.zeros(ps)
    for k in range(ps):
        x = np.zeros(n)
        for j in range(n):
            bits = Pn[k, nb*j : nb*(j+1)]
            adec = int("".join(map(str, bits)), 2)
            x[j] = xl[j] + (xu[j] - xl[j]) / (2**nb - 1) * adec
        ft[k] = fun(x)
    return ft

def frank(ps, ft):
    indices = np.arange(1, ps + 1)
    # ft 기준으로 내림차순 정렬
    sorted_idx = np.argsort(ft)[::-1]
    ft_sorted = ft[sorted_idx]
    frk = indices[sorted_idx]
    return frk, ft_sorted

def roul(ps, am):
    temp = np.random.rand()
    for k in range(ps):
        if temp <= am[k]:
            return k
    return ps - 1

def cros(Pn, nb, n, kstr1, kstr2):
    kchd = np.zeros(nb * n, dtype=int)
    for j in range(n):
        jp = nb * j
        nc = np.random.randint(0, nb + 1)
        if nc == 0:
            kchd[jp:jp+nb] = Pn[kstr2, jp:jp+nb]
        elif nc == nb:
            kchd[jp:jp+nb] = Pn[kstr1, jp:jp+nb]
        else:
            kchd[jp:jp+nc] = Pn[kstr1, jp:jp+nc]
            kchd[jp+nc:jp+nb] = Pn[kstr2, jp+nc:jp+nb]
    return kchd

def mutat(kchd, mp, nb, n):
    mask = np.random.rand(nb * n) <= mp
    kchd[mask] = 1 - kchd[mask]
    return kchd