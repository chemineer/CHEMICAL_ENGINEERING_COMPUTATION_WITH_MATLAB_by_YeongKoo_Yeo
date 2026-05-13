import numpy as np
from scipy.interpolate import interp1d, PchipInterpolator, CubicSpline

# 1. 데이터 정의 (B-T 평형 데이터)
x = np.array([0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0])
y = np.array([0.00, 0.21, 0.38, 0.51, 0.62, 0.71, 0.79, 0.86, 0.91, 0.96, 1.00])

# 보간을 수행할 목표 지점
xi = 0.45

# 2. 다양한 보간 방법 수행 (MATLAB interp1 대응)

# 선형 보간 (Linear interpolation)
# MATLAB: interp1(x, y, 0.45, 'linear')
f_linear = interp1d(x, y, kind='linear')
lv = f_linear(xi)

# 조각적 삼차 에르미트 보간 (Piecewise cubic Hermite interpolation)
# MATLAB: interp1(x, y, 0.45, 'pchip')
f_pchip = PchipInterpolator(x, y)
pv = f_pchip(xi)

# 조각적 삼차 스플라인 보간 (Piecewise cubic spline)
# MATLAB: interp1(x, y, 0.45, 'spline')
f_spline = CubicSpline(x, y) # 또는 interp1d(x, y, kind='cubic')
sv = f_spline(xi)

# 최근접 이웃 보간 (Nearest neighbor interpolation)
# MATLAB: interp1(x, y, 0.45, 'nearest')
f_nearest = interp1d(x, y, kind='nearest')
nv = f_nearest(xi)

# 3. 결과 출력
print(f"linear:  {lv:6.4f}")
print(f"pchip:   {pv:6.4f}")
print(f"spline:  {sv:6.4f}")
print(f"nearest: {nv:6.4f}")