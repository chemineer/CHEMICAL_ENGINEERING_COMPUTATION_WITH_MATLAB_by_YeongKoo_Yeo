import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

# 1. 데이터 정의 (실험 데이터)
# x: Flow rate (liter/sec), y: Pressure drop (kPa)
x = np.array([0, 9.7, 15.6, 21.3, 31.7, 35.2, 38.4, 42.9])
y = np.array([0, 0.28, 0.524, 0.998, 1.695, 2.306, 2.781, 3.205])

# 2. 삼차 스플라인 보간 수행
# MATLAB의 spline(x, y, xi)와 동일한 기능을 위해 CubicSpline 객체 생성
cs = CubicSpline(x, y)

# 3. 보간할 지점 설정 (min(x)부터 max(x)까지 0.1 간격)
xi = np.arange(x.min(), x.max() + 0.1, 0.1)
yi = cs(xi)

# 4. 시각화
plt.figure(figsize=(8, 5))
plt.plot(xi, yi, label='Cubic spline interpolation') # 보간된 곡선
plt.plot(x, y, 'o', label='Experimental data')      # 원본 데이터 포인트

plt.xlabel('Flow rate(liter/sec)')
plt.ylabel('Pressure drop(kPa)')
plt.legend(loc='best')
plt.grid(True)
plt.title('Cubic Spline Interpolation of Pressure Drop Data')
plt.show()