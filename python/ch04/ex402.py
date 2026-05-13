import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# 상수 정의
B = -0.0361
C = 2.7047e-3
D = -4.4944e-4
R = 0.08206
T = 200
n = 300

# 압력 범위 생성 (1부터 30 atm까지 300개 지점)
P = np.linspace(1, 30, n)

# 초기 추정값 (x0)
x0 = 0.5 * np.ones(n)

# 비선형 방정식 정의 (x는 밀도 rho를 의미)
# 식: 1 + B*x + C*x^2 + D*x^3 - P/(x*R*T) = 0
def f(x, p_val):
    return 1 + B * x + C * x**2 + D * x**3 - p_val / (x * R * T)

# 각 압력(P)에 대해 fsolve를 사용하여 밀도(rho) 계산
rho = np.zeros(n)
for i in range(n):
    # 각 P[i]에 대응하는 방정식의 해를 구함
    rho[i] = fsolve(f, x0[i], args=(P[i],))[0]

# 결과 그래프 그리기
plt.figure(figsize=(8, 5))
plt.plot(P, rho)
plt.xlabel('P(atm)')
plt.ylabel('Density(mol/liter)')
plt.title('Density vs Pressure')
plt.grid(True)
plt.show()