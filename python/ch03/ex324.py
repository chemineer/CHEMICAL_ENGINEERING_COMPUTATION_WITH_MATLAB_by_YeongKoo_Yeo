import numpy as np
import matplotlib.pyplot as plt
from hvapn import hvapn

# 온도 범위 설정 (0부터 350까지)
T = np.arange(0, 351)

# 기화열 계산
hv = hvapn(T, 'H2O')

# 그래프 그리기
plt.figure(figsize=(10, 6))
plt.plot(T, hv)
plt.xlabel('T(C)')
plt.ylabel('Heat of vaporization of water(cal/g)')
plt.title('Heat of Vaporization of Water')
plt.grid(True)
plt.show()