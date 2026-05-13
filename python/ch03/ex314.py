import numpy as np
import matplotlib.pyplot as plt
from condG import condG

# 1. 온도 범위 설정 (25도에서 900도까지, 켈빈 온도로 변환)
# MATLAB: T = [25:900] + 273.15
t_celsius = np.arange(25, 901)
T = t_celsius + 273.15

# 2. 기체 열전도도 계산 함수 호출
# MATLAB: k = condG(T, 'propane')
k = condG(T, 'propane')

# 3. 결과 시각화 (MATLAB의 plot, xlabel, axis tight, grid 대응)
plt.figure(figsize=(8, 6))
plt.plot(t_celsius, k)
plt.xlabel('T(C)')
plt.ylabel('Thermal conductivity of propane(μcal/s/cm/K)')
plt.title('Thermal Conductivity of Propane vs Temperature')

# MATLAB의 axis tight 대응 (데이터 범위에 맞춰 축 축소)
plt.autoscale(enable=True, axis='x', tight=True)
plt.grid(True)

plt.show()