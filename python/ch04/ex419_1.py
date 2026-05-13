import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# 1. 데이터 정의
A1, B1, C1 = 14.2724, 2945.47, 224.0
A2, B2, C2 = 14.2043, 2972.64, 209.0
P = 70.0  # 총 압력

# Antoine 식을 이용한 포화압력 함수
def psat(A, B, C, T):
    return np.exp(A - B / (T + C))

# 2. 루프 계산
x1_range = np.linspace(0, 1, 101)
T_res = []
y1_res = []

T0 = 60.0 # 초기 온도 추정값

for x1 in x1_range:
    # 방정식 정의: f(T) = x1*Psat1(T) + (1-x1)*Psat2(T) - P = 0
    # MATLAB 코드의 fp 식을 라울의 법칙 기준으로 재구성
    def func(T):
        return x1 * psat(A1, B1, C1, T) + (1 - x1) * psat(A2, B2, C2, T) - P
    
    # fsolve로 온도 T 구하기
    T_sol = fsolve(func, T0)
    T_res.append(T_sol[0])
    
    # 기상 조성 y1 계산: y1 = (x1 * Psat1) / P
    y1 = (x1 * psat(A1, B1, C1, T_sol[0])) / P
    y1_res.append(y1)
    
    # 다음 반복을 위해 T0 업데이트
    T0 = T_sol[0]

# 3. 시각화
plt.figure(figsize=(8, 6))
plt.plot(x1_range, T_res, label='T-x1 (Bubble Point Curve)')
plt.plot(y1_res, T_res, '.-', label='T-y1 (Dew Point Curve)')
plt.xlabel('x1, y1')
plt.ylabel('T (°C)')
plt.title(f'T-xy Diagram at P={P} kPa')
plt.legend(loc='best')
plt.grid(True)
plt.show()