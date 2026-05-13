import numpy as np

def jacob(fun, x, h):
    """
    야코비 행렬(Jacobian)의 수치적 근사 계산
    :param fun: 함수 (입력 x에 대해 결과값을 반환하는 함수)
    :param x: 현재 지점 (리스트 또는 numpy 배열)
    :param h: 미분 간격 (step size)
    :return: Hs (야코비 행렬)
    """
    x = np.array(x, dtype=float).flatten()
    n = len(x)
    hd = 2 * h
    M = np.eye(n)
    
    # 결과 행렬 초기화 (출력 크기 확인을 위해 먼저 함수 한번 호출)
    f_sample = np.array(fun(x))
    m_out = len(f_sample)
    Hs = np.zeros((m_out, n))
    
    # 중앙 차분법(Central Difference)을 이용한 야코비 행렬 계산
    for k in range(n):
        # x + h*ek 와 x - h*ek 계산
        x_plus = x + M[k, :] * h
        x_minus = x - M[k, :] * h
        
        # 각 방향에 대한 미분 근사
        # (f(x+h) - f(x-h)) / 2h
        Hs[:, k] = (np.array(fun(x_plus)) - np.array(fun(x_minus))) / hd
        
    return Hs