import numpy as np
import matplotlib.pyplot as plt

# 1. 원래 함수 f(x) 정의
f = lambda x: 0.3 + 20*x - 180*x**2 + 650*x**3 - 880*x**4 + 360*x**5

# 2. 수치 미분을 위한 x 범위 설정 (간격 0.1)
x = np.arange(0, 1.1, 0.1)
y = f(x)

# 3. 수치 미분 수행 (dy/dx)
# numpy.diff(y)는 y[i+1] - y[i]를 계산합니다.
dr = np.diff(y) / np.diff(x)

# 4. 구간의 중점(Midpoint) 계산
# diff 결과는 원래 데이터보다 크기가 1 작으므로, 비교를 위해 x의 중점 값을 구합니다.
xm = (x[:-1] + x[1:]) / 2

# 5. 실제 도함수(Exact differentiation) 정의 및 계산
# f'(x) = 20 - 360*x + 1950*x^2 - 3520*x^3 + 1800*x^4
xp = np.arange(0, 1.01, 0.01)
yp = 20 - 360*xp + 1950*xp**2 - 3520*xp**3 + 1800*xp**4

# 6. 결과 시각화
plt.figure(figsize=(8, 5))
plt.plot(xp, yp, label='Exact differentiation', color='blue')
plt.plot(xm, dr, 'ro', label='Numerical differentiation') # 'ro'는 빨간색 원형 마커

plt.xlabel('x')
plt.ylabel('df(x)/dx')
plt.title('Comparison of Exact and Numerical Differentiation')
plt.legend(loc='best')
plt.grid(True)
plt.show()