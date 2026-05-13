import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from adfun import adfun

# 매개변수 설정
P = 162
Vspan = [0, 4]
FN2 = 28.3
X0 = [38.3 - FN2, 0, 0, 1150]
pf = (P, FN2)

# solve_ivp를 이용한 적분 (MATLAB의 ode45 대응)
# t_eval을 통해 원하는 출력 지점을 지정합니다.
sol = solve_ivp(lambda V, X: adfun(V, X, pf), Vspan, X0, method='RK45', t_eval=np.linspace(Vspan[0], Vspan[1], 100))

V = sol.t
X = sol.y.T  # 행과 열을 맞추기 위해 전치

# 결과 계산
xc = (X0[0] - X[:, 0]) / X0[0]

# 그래프 출력
plt.figure(figsize=(12, 4))

plt.subplot(1, 2, 1)
plt.plot(V, X[:, 3]) # MATLAB X(:,4)는 파이썬 인덱스 3
plt.xlabel('Reactor volume(m^3)')
plt.ylabel('Temperature(K)')
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(V, xc)
plt.xlabel('Reactor volume(m^3)')
plt.ylabel('Conversion')
plt.grid(True)

plt.tight_layout()
plt.show()