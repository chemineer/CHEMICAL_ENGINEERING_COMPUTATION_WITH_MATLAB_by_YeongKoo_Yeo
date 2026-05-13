import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from pfrconc import pfrconc

# 1. 데이터 및 파라미터 설정 (MATLAB의 pf 구조체 대응)
n = 20
pf = {
    'k': 0.18,
    'v': 0.5,
    'C0': 1,
    'L': 0.5,
    'n': n
}
h = pf['L'] / n
C0 = pf['C0']

# 2. 초기 조건 및 시간 범위 설정
Z0 = np.ones(n) * C0
tspan = (0, 10)
t_eval = np.linspace(0, 10, 100) # 부드러운 그래프를 위한 시간 타임스탬프

# 3. ODE solver (ode45 대신 solve_ivp 사용)
# solve_ivp는 함수 인자 순서가 (t, y)이며, 추가 인자는 args로 전달합니다.
sol = solve_ivp(pfrconc, tspan, Z0, args=(pf,), t_eval=t_eval)

t = sol.t
C = sol.y.T # (시간, 위치) 형태로 전치

# 4. 결과 처리
# MoL에 의한 정상 상태 (마지막 시간대의 농도)
Cs = np.insert(C[-1, :], 0, C0)

# 정확한 해(Exact Solution) 계산
x = np.arange(0, pf['L'] + h, h)
Cm = C0 * np.exp(-pf['k'] * x / pf['v'])

# 5. 시각화 (MATLAB의 subplot 구조 재현)
plt.figure(figsize=(12, 5))

# Subplot 1: 정상 상태 비교 (MoL vs Exact)
plt.subplot(1, 2, 1)
plt.plot(x, Cs, 'b-', label='St-st by MoL')
plt.plot(x, Cm, 'r--', label='St-st by exact solution')
plt.xlabel('x(m)')
plt.ylabel('C(mol/liter)')
plt.legend()
plt.grid(True)

# Subplot 2: 출구 농도(C_end)의 시간 변화
plt.subplot(1, 2, 2)
plt.plot(t, C[:, -1], 'g-')
plt.xlabel('t(min)')
plt.ylabel('C(mol/liter)')
plt.title('Concentration at Exit vs Time')
plt.grid(True)

plt.tight_layout()
plt.show()