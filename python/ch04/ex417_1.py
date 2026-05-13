import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from binflash import binflash

# 1. 환경 설정 및 데이터 정의
A = np.array([13.8183, 13.8587])
B = np.array([2477.07, 2991.32])
C = np.array([233.21, 216.64])
T = 60.0
z = 0.6
v_range = np.linspace(0, 1, 101)  # 0부터 1까지 0.01 간격으로 101개 포인트
P_list = []

# 2. 기화율 v에 따른 평형 압력 P 계산
x0 = [0.1, 0.6, 50.0]  # 초기 추정값 (x1, y1, P)

for v in v_range:
    # fsolve로 방정식 해 구하기
    sol = fsolve(binflash, x0, args=(v, T, z, A, B, C))
    P_list.append(sol[2])  # 결과 중 압력(P) 값 저장
    x0 = sol  # 다음 반복을 위해 현재 해를 초기값으로 사용 (수렴 속도 향상)

# 3. 결과 시각화
plt.figure(figsize=(8, 5))
plt.plot(v_range, P_list, label='Equilibrium Pressure')
plt.xlabel('Vaporized fraction (v)')
plt.ylabel('Pressure (kPa)')
plt.title('Vaporized Fraction vs. Equilibrium Pressure')
plt.grid(True)
plt.legend()
plt.show()