import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# 상수 및 물성치 설정
D = 0.05            # 지름 (m)
h = 98.6            # 대류 열전달 계수 (W/m^2·K)
k = 1.18            # 열전도도 (W/m·K)
alpa = 4.97e-7      # 열확산율 (m^2/s)
Ta = 25             # 주변 온도 (C)
Ti = 340            # 초기 온도 (C)

# 시간 배열 설정 (0분에서 15분까지 2000개의 지점)
t = np.linspace(0, 15 * 60, 2000) # 초 단위 변환

# 초월 방정식 정의: x*tan(x) - Bi = 0
# Bi (Biot number) 관련 항: h * (D/2) / k
f = lambda x: x * np.tan(x) - h * D / k / 2 #

# 방정식의 해(gamma) 찾기 (초기 추정치 x0 = 0.1)
x0 = 0.1 #
gam = fsolve(f, x0)[0] #

# 온도 계산 식 적용
# T = Ta + (Ti-Ta) * [4*sin(gam) * exp(-4*gam^2*alpa*t/D^2) / (2*gam + sin(2*gam))]
term_top = 4 * np.sin(gam) * np.exp(-4 * gam**2 * alpa * t / D**2) #
term_bottom = 2 * gam + np.sin(2 * gam) #
T = Ta + (Ti - Ta) * term_top / term_bottom #

# 결과 출력
print(f"gamma = {gam:g}") #

# 그래프 시각화
plt.figure(figsize=(8, 5)) #
plt.plot(t / 60, T) # x축을 분(min) 단위로 변환하여 출력
plt.grid(True) #
plt.xlabel('t(min)') #
plt.ylabel('T(deg.C)') #
plt.title('Unsteady-state Heat Conduction') #
plt.show() #