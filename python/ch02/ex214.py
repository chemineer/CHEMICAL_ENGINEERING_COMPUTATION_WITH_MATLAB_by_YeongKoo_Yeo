import numpy as np
from scipy.optimize import fsolve

# 1. 안토완 상수 (A, B, C) 설정
# 인덱스 0: Benzene, 1: Toluene
A = np.array([15.90085, 16.01066])
B = np.array([2788.507, 3094.543])
C = np.array([220.790, 219.377])

# 2. 고정 조건 설정
P_total = 760  # 전체 압력 (mmHg)
x = np.array([0.4, 0.6])  # 액상 몰 분율 (x1, x2)

# 3. 버블 포인트 온도 계산을 위한 함수 정의
# f(t) = Σ(xi * Pi_sat) - P_total = 0
def bubble_point_func(t):
    # 각 성분의 증기압 계산 (안토완 방정식)
    # P_sat = exp(A - B / (t + C))
    p_sat1 = np.exp(A[0] - B[0] / (t + C[0]))
    p_sat2 = np.exp(A[1] - B[1] / (t + C[1]))
    
    # 혼합물 전체 압력과 설정 압력의 차이 반환
    return x[0] * p_sat1 + x[1] * p_sat2 - P_total

# 4. 수치해석을 통한 온도(t) 도출
t_guess = 50  # 초기 추측값 (50 deg.C)
t_solution = fsolve(bubble_point_func, t_guess)

print(f"The bubble point temperature = {t_solution[0]:.4f} deg.C")