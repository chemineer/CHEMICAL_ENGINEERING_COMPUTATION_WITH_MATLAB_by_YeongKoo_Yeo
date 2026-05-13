import numpy as np
from scipy.interpolate import RectBivariateSpline

# 1. 데이터 정의 (MATLAB과 동일)
# x: [2, 10], y: [1, 8]
x = np.array([2, 10])
y = np.array([1, 8])

# z 데이터: MATLAB의 [80 78; 75 90]은 행(y), 열(x) 구조입니다.
# MATLAB z(1,1)=80, z(1,2)=78, z(2,1)=75, z(2,2)=90
z = np.array([[80, 78], 
              [75, 90]])

# 2. 2차원 Spline 보간 객체 생성
# kx, ky=1은 데이터 포인트가 적을 때 선형으로 작동하며, 
# 데이터가 많을 경우 kx, ky=3을 사용하여 Cubic Spline을 구현합니다.
interp_func = RectBivariateSpline(x, y, z.T, kx=1, ky=1)

# 3. 특정 지점(xi=6.4, yi=5.2)에서의 값 계산
xi = 6.4
yi = 5.2
zi = interp_func(xi, yi)

print(f"보간된 결과값 (zi): {zi[0][0]:.4f}")