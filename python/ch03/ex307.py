import numpy as np
import matplotlib.pyplot as plt
from gunyam import gunyam

# 1. 벤젠(Benzene) 물성치 설정
Tc = 288.93   # 임계 온도 (deg.C)
Pc = 49.24    # 임계 압력 (bar)
w = 0.212     # 편심 인자
Mw = 78       # 분자량 (g/mol)

# 2. 온도 범위 설정 (MATLAB: T = 0:70)
T = np.arange(0, 71)  # 0부터 70까지의 정수 배열
n = len(T)

# 3. 루프를 통한 밀도 계산
denBz = np.zeros(n)
for k in range(n):
    # gunyam 함수를 사용하여 밀도 계산
    denBz[k] = gunyam(Tc, Pc, w, Mw, T[k])

# 4. 시각화 (MATLAB의 plot, grid, xlabel, ylabel 대응)
plt.figure(figsize=(8, 5))
plt.plot(T, denBz)
plt.grid(True)
plt.xlabel('T(deg.C)')
plt.ylabel('Density(kg/m^3)')
plt.title('Density of liquid benzene')
plt.show()