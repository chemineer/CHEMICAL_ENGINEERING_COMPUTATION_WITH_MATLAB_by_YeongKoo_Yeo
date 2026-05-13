import numpy as np

# 1. 함수 및 도함수 정의 (람다 함수 이용)
f = lambda x: 0.3 + 20*x - 180*x**2 + 650*x**3 - 880*x**4 + 360*x**5
df = lambda x: 20 - 360*x + 1950*x**2 - 3520*x**3 + 1800*x**4

# 2. x 데이터 생성 (0부터 1.0까지 0.1 간격)
x = np.arange(0, 1.1, 0.1)
y = f(x)

# 3. 도함수 계산
# np.gradient(y, spacing): y값들 사이의 간격(0.1)을 인자로 전달합니다.
dr = np.gradient(y, 0.1)  # 수치적 구배(Estimation by gradient)
dy = df(x)               # 정확한 해(Exact solution)

# 4. 결과 출력
print("x 값:", x)
print("-" * 60)
print(f"{'x':>5} | {'dr (Numerical)':>15} | {'dy (Exact)':>15} | {'Error':>10}")
print("-" * 60)

for i in range(len(x)):
    error = abs(dr[i] - dy[i])
    print(f"{x[i]:5.1f} | {dr[i]:15.6f} | {dy[i]:15.6f} | {error:10.6f}")