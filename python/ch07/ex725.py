import numpy as np
from scipy.integrate import quad

# 데이터 정의
x = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
y = np.array([0.000, 0.211, 0.378, 0.512, 0.623, 0.714, 0.791, 0.856, 0.911, 0.959, 1.000])

# 4차 다항식 피팅 (MATLAB의 polyfit과 동일)
p = np.polyfit(x, y, 4)

# 피팅된 다항식 함수 정의 (p[0]*x^4 + p[1]*x^3 + ...)
poly_func = np.poly1d(p)

# 적분 대상 함수 f(x) 정의
def f(x_val):
    return 1.0 / (x_val - poly_func(x_val))

# 초기 설정값
x1 = 0.5833
vfix = 0.8755
crit = 1e-2
xf = None

# x2를 x1부터 0.0까지 -1e-3 간격으로 반복 (MATLAB의 for 루프 재현)
# np.arange는 끝점(0.0)을 포함하기 위해 약간 더 작은 값을 설정합니다.
for x2 in np.arange(x1, -0.001, -0.001):
    # scipy의 quad 함수를 이용한 수치 적분 (MATLAB의 integral과 동일)
    res, err = quad(f, x1, x2)
    
    if abs(res - vfix) <= crit:
        xf = x2
        break

print(f"xf = {xf:.4f}")