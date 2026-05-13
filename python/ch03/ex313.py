import numpy as np
import matplotlib.pyplot as plt
from condL import condL

# 1. 온도 범위 설정 (0도에서 350도까지, 켈빈 온도로 변환)
# MATLAB: T = [0:350] + 273.15
t_celsius = np.arange(0, 351)
T = t_celsius + 273.15

# 2. 열전도도 계산 함수 호출
# MATLAB: k = condL(T, 'water')
k = condL(T, 'water')

# 3. 결과 시각화 (MATLAB의 plot, xlabel, ylabel, grid 대응)
plt.figure(figsize=(8, 6))
plt.plot(t_celsius, k)
plt.xlabel('T(C)')
plt.ylabel('Thermal conductivity of water(μcal/s/cm/C)')
plt.title('Thermal Conductivity of Water vs Temperature')
plt.grid(True)
plt.show()