import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import fsolve

# 1. 데이터 및 초기 파라미터 설정
x1i, x2i, x3i = 1.2e-3, 2e-3, 2.5e-3
Kd = 1e-7
Ke = 1e-14
Fa = 0.01667
Fb = 2.333e-3
V = 2.5
av = Fa / V
bv = Fb / V
xi = [x1i, x2i, x3i]

# 2. ODE 시스템 정의
def dx_dt(x, t):
    x1, x2, x3 = x
    dx1 = av * (x1i - x1) - bv * x1
    dx2 = bv * (x2i - x2) - av * x2
    dx3 = bv * (x3i - x3) - av * x3
    return [dx1, dx2, dx3]

# 3. ODE 풀기 (x축 0~250 범위 설정)
t = np.linspace(0, 250, 500)
sol = odeint(dx_dt, xi, t)

x1 = sol[:, 0]
x2 = sol[:, 1]
x3 = sol[:, 2]

# 4. pH 계산 ([H+] 구하기)
# fsolve를 사용하여 각 시간대별 비선형 방정식의 해를 구함
def calculate_pH(x1_val, x2_val, x3_val):
    # h = [H+]
    f = lambda h: h + x2_val + x3_val - x1_val - Ke/h - x3_val/(1 + (Kd * h / Ke))
    
    # 초기 추정값: 중성에 가까운 값 혹은 이전 pH 기반 (안정성을 위해 1e-7 근처 설정)
    hp_initial = 1e-7
    hp_solution = fsolve(f, hp_initial)
    return -np.log10(np.abs(hp_solution[0]))

# 각 타임스텝에 대해 pH 계산
pH_values = [calculate_pH(x1[i], x2[i], x3[i]) for i in range(len(t))]

# 5. 시각화
plt.figure(figsize=(12, 5))

# 농도 그래프
plt.subplot(1, 2, 1)
plt.plot(t, x1, label='$x_1$', color='blue')
plt.plot(t, x2, label='$x_2$', linestyle=':', color='green')
plt.plot(t, x3, label='$x_3$', linestyle='--', color='red')
plt.grid(True, linestyle='--', alpha=0.7)
plt.xlabel('t (sec)')
plt.ylabel('$x_i$ (mol/liter)')
plt.title('Concentrations over Time')
plt.legend()

# pH 그래프
plt.subplot(1, 2, 2)
plt.plot(t, pH_values, color='magenta')
plt.grid(True, linestyle='--', alpha=0.7)
plt.xlabel('t (sec)')
plt.ylabel('pH')
plt.title('pH Profile (0-250s)')

plt.tight_layout()
plt.show()