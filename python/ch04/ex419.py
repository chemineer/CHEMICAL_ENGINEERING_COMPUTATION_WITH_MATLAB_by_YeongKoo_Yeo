import numpy as np
import matplotlib.pyplot as plt

# 1. 데이터 정의
T = 75.0
# Antoine 상수 (성분 1, 2)
A1, B1, C1 = 14.2724, 2945.47, 224.0
A2, B2, C2 = 14.2043, 2972.64, 209.0

# Antoine 식을 이용한 포화압력 계산 함수
def psat(A, B, C, T):
    return np.exp(A - B / (T + C))

Psat1 = psat(A1, B1, C1, T)
Psat2 = psat(A2, B2, C2, T)

# 2. 조성 x1에 따른 P 및 y1 계산
x1 = np.linspace(0, 1, 101)  # 0부터 1까지 101개 포인트
# 라울의 법칙: P = x1*Psat1 + (1-x1)*Psat2
P = x1 * Psat1 + (1 - x1) * Psat2
# 기상 조성 y1: y1 = (x1 * Psat1) / P
y1 = (x1 * Psat1) / P

# 3. 시각화
plt.figure(figsize=(8, 6))
plt.plot(x1, P, label='P-x1')
plt.plot(y1, P, '.-', label='P-y1')
plt.xlabel('x1, y1')
plt.ylabel('P')
plt.title(f'P-xy Diagram at T={T}°C')
plt.legend(loc='best')
plt.grid(True)
plt.show()