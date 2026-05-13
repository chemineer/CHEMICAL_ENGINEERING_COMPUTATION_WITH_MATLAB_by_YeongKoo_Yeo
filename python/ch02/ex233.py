import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

# 1. 데이터 정의
tm = np.array([0, 18, 42, 55, 70, 82, 86, 95, 102, 115]) # 시간(min)
TC = np.array([15, 23, 24, 36, 78, 80, 98, 96, 127, 126]) # 온도(deg.C)

# 2. PCHIP 보간 수행
# MATLAB: yp = interp1(tm, TC, tinv, 'pchip')
tinv = np.linspace(0, 115, 100) # 보간할 시간 범위 생성
pchip_func = PchipInterpolator(tm, TC) # PCHIP 보간 객체 생성
yp = pchip_func(tinv) # 보간값 계산

# 3. 시각화
plt.figure(figsize=(9, 6))
plt.plot(tm, TC, 'o', label='Data') # 원본 데이터 포인트
plt.plot(tinv, yp, label='Cubic Hermite interpolation') # 보간 곡선

plt.xlabel('Time(min)')
plt.ylabel('Temperature(deg.C)')
plt.axis([0, 115, 0, 130]) # 축 범위 설정
plt.title('1D Interpolation using PCHIP Method')
plt.legend(loc='best')
plt.grid(True)
plt.show()