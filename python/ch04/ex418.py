import numpy as np
from scipy.optimize import fsolve

# 1. 데이터 정의
# Antoine 상수 (로그 밑이 10인 경우: log10(P_sat) = A - B / (C + T_C))
A = np.array([4.00139, 4.02023, 4.05075, 4.07356])
B = np.array([1170.875, 1263.909, 1356.360, 1438.03])
C = np.array([224.317, 216.432, 209.635, 202.694])

x = np.array([0.32, 0.31, 0.25, 0.12])  # 성분 조성
P = 1.5  # 전체 압력 (bar)

# 2. 비선형 방정식 정의
# 라울의 법칙: sum(x_i * P_sat_i) = P_total
# 식: sum(x_i * 10^(A - B / (C + T_C))) / P - 1 = 0
def bubble_point_eq(T_C):
    P_sat = 10**(A - B / (C + T_C))
    return np.sum(x * P_sat) / P - 1

# 3. 수치 해석 (fsolve)
T0 = 400 - 273.15  # 초기 추정값 (섭씨 온도 기준 변환)
T_solution = fsolve(bubble_point_eq, T0)

print(f"버블 포인트 온도: {T_solution[0]:.2f} °C")