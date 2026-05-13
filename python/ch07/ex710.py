import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 외부 모듈 임포트
from vapdat import * # 데이터 정의 (변수들이 포함되어 있어야 함)
from vapr import vapr    # 미분 방정식 함수 vapr(t, z) 정의

# 1. 초기 조건 및 시간 범위 설정
# MATLAB: mL0 = 2800; mV0 = 100; z0 = [mV0 mL0];
mL0 = 2800
mV0 = 100
z0 = [mV0, mL0]
t_span = (0, 0.1)

# 2. ODE 풀이 (ode45 -> solve_ivp)
# vapr.vapr는 vapr.py 파일 내에 def vapr(t, z): 형식으로 정의되어 있어야 합니다.
sol = solve_ivp(
    vapr, 
    t_span, 
    z0, 
    method='RK45', 
    t_eval=np.linspace(t_span[0], t_span[1], 100)
)

# 3. 결과 시각화
plt.figure(figsize=(8, 5))
plt.plot(sol.t, sol.y[0, :], label='$m_V$')       # z(:,1) -> sol.y[0]
plt.plot(sol.t, sol.y[1, :], ':', label='$m_L$')  # z(:,2) -> sol.y[1]

plt.grid(True)
plt.xlabel('t(h)')
plt.ylabel('$m_V, m_L(kg)$')
plt.legend()
plt.show()