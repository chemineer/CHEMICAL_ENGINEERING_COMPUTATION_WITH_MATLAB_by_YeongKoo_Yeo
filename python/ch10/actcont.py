import numpy as np
from typing import Callable

def actcont(n: int, m: int, ne: int, crit: float, x: np.ndarray, A: np.ndarray, b: np.ndarray):
    """
    활동 제약 조건(Active Constraints)을 관리하는 함수
    """
    # g = A*x' - b 형태를 맞추기 위해 x를 열 벡터로 변환
    g = A @ x.T - b
    
    # nca 초기화 (1부터 m까지)
    nca = np.arange(1, m + 1)
    nc = ne
    
    # 제약 조건 검사
    for j in range(ne, m):
        if abs(g[j]) < crit:
            nc += 1
            # nca 배열 요소 교환
            ntemp = nca[j]
            nca[j] = nca[nc - 1] # 파이썬 인덱스는 0부터 시작하므로 nc-1
            nca[nc - 1] = ntemp
            
    return nc, nca

def bisec(delfun: Callable[[np.ndarray], np.ndarray], n: int, alphak: float, x0: np.ndarray, d: np.ndarray):
    """
    이분법(Bisection method)을 이용한 선 탐색(Line search) 함수
    """
    mcrit = 1e-6
    a1 = 0.0
    a2 = alphak
    aw = a2 - a1
    
    while (a2 - a1) > mcrit * aw:
        am = (a1 + a2) / 2.0
        x_curr = x0 + am * d
        
        # delfun(x)' * d' 계산
        # 파이썬에서는 배열 연산이므로 내적(dot product)으로 처리
        fp = np.dot(delfun(x_curr), d)
        
        if fp < 0:
            a1 = am
        elif fp > 0:
            a2 = am
        else:
            break
            
    x = x0 + a1 * d
    return x