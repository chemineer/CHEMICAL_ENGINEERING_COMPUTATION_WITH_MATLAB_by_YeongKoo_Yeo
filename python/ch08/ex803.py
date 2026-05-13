import numpy as np
from scipy.optimize import fsolve

def frad(x, radt):
    """
    복사 열전달 평형 방정식을 정의하는 함수
    x[0] = J1, x[1] = J2, x[2] = T2
    """
    # 구조체 데이터 추출
    A1, A2 = radt['A1'], radt['A2']
    W, H, L = radt['W'], radt['H'], radt['L']
    sig, T1, T3 = radt['sig'], radt['T1'], radt['T3']
    ep1, ep2, q2_const = radt['ep1'], radt['ep2'], radt['q2']
    
    J1, J2, T2 = x[0], x[1], x[2]
    
    # 형태 계수(View Factor) 계산
    X_dim = W / H
    Y_dim = L / H
    
    # F12 계산 (복잡한 로그 및 삼각함수 항 포함)
    term1 = np.log(np.sqrt((1 + X_dim**2) * (1 + Y_dim**2) / (1 + X_dim**2 + Y_dim**2)))
    term2 = X_dim * np.sqrt(1 + Y_dim**2) * np.arctan(X_dim / np.sqrt(1 + Y_dim**2))
    term3 = Y_dim * np.sqrt(1 + X_dim**2) * np.arctan(Y_dim / np.sqrt(1 + X_dim**2))
    term4 = X_dim * np.arctan(X_dim)
    term5 = Y_dim * np.arctan(Y_dim)
    
    F12 = (2 / (np.pi * X_dim * Y_dim)) * (term1 + term2 + term3 - term4 - term5)
    F13 = 1 - F12
    F23 = F13 * A1 / A2
    
    # 비선형 방정식 시스템 정의
    f1 = (ep1 * A1 * (sig * T1**4 - J1) / (1 - ep1)) - (F12 * A1 * (J1 - J2)) - (F13 * A1 * (J1 - sig * T3**4))
    f2 = (ep2 * A2 * (sig * T2**4 - J2) / (1 - ep2)) - (F12 * A1 * (J2 - J1)) - (F23 * A2 * (J2 - sig * T3**4))
    f3 = q2_const + (ep2 * A2 * (sig * T2**4 - J2) / (1 - ep2))
    
    return [f1, f2, f3]

# 1. 파라미터 설정 (딕셔너리 사용)
radt = {
    'A1': 10, 'A2': 15, 'W': 1, 'H': 1, 'L': 10,
    'sig': 5.67e-8, 'T1': 1000, 'T3': 300, 'ep1': 0.9,
    'ep2': 0.5, 'q2': 77100
}

# 2. 초기 추정값 설정 및 방정식 풀이
x0 = [500, 500, 500]
x_solution = fsolve(frad, x0, args=(radt,))

# 3. 결과 추출 및 추가 열유속(q) 계산
J1_sol, J2_sol, T2_sol = x_solution

q1 = radt['ep1'] * radt['A1'] * (radt['sig'] * (radt['T1']**4) - J1_sol) / (1 - radt['ep1'])
q2_check = -radt['ep2'] * radt['A2'] * (radt['sig'] * T2_sol**4 - J2_sol) / (1 - radt['ep2'])

# 4. 결과 출력
print(f"J1 = {J1_sol:.4f}, J2 = {J2_sol:.4f}, T2 = {T2_sol:.4f}")
print(f"q1 = {q1:.4f}")
print(f"q2 (calculated) = {q2_check:.4f}")