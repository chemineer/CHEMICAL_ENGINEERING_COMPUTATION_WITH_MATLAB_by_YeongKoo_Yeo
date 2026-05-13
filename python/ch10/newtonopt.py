import numpy as np
from jacob import jacob

def newtonopt(fun, x0, crit, kmax):
    """
    뉴턴법(Newton's method)을 이용한 최적화
    :param fun: 목적 함수
    :param x0: 시작점 (리스트 또는 numpy 배열)
    :param crit: 정지 기준 (허용 오차)
    :param kmax: 최대 반복 횟수
    :return: xopt, fopt, iter
    """
    h = 1e-4
    x = np.array(x0, dtype=float).flatten()
    fx = np.array(fun(x))
    
    # 반복문 시작
    for k in range(kmax):
        # 야코비 행렬(jacob)을 사용하여 dx 계산
        # dx = -J^-1 * f(x)
        # 파이썬에서는 numpy.linalg.solve를 사용하는 것이 더 안정적입니다.
        J = jacob(fun, x, h)
        dx = -np.linalg.solve(J, fx)
        
        # 해 갱신
        x = x + dx
        fx = np.array(fun(x))
        
        # 정지 조건 검사
        if np.linalg.norm(fx) < crit or np.linalg.norm(dx) < crit:
            return x, fx, k + 1
            
    return x, fx, kmax

# 참고: 위 코드 실행을 위해서는 이전에 변환한 jacob 함수가 필요합니다.