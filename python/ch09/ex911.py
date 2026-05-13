import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from LTmodel import LTmodel  # LTmodel.py 파일에서 함수 임포트

# 1. 데이터 구조체 설정 (매트랩의 ht 구조체 대응)
class Struct: pass
ht = Struct()
ht.A1 = 0.25; ht.A2 = 0.25; ht.F0 = 0.4; ht.T0 = 25; ht.H = 0.5
ht.c1 = 0.6; ht.c2 = 0.6; ht.rCp = 4180; ht.Q1 = 6000; ht.Q2 = 6000

# 2. 초기 조건 및 시뮬레이션 시간 설정
h10 = 0.4; h20 = 0.35
z0 = [h10, h20, 25, 25]  # [h1, h2, T1, T2]
t_span = (0, 10)         # tspan = [0 10]
t_eval = np.linspace(0, 10, 500) # 그래프를 위한 시간 간격

# 3. ODE 풀이 수행 (ode45 대응)
# LTmodel 함수에 추가 인자 ht를 전달하기 위해 args를 사용합니다.
sol = solve_ivp(
    fun=LTmodel,
    t_span=t_span,
    y0=z0,
    args=(ht,),
    t_eval=t_eval,
    method='RK45'
)

# 결과 추출
t = sol.t
h1, h2, T1, T2 = sol.y

# 4. 결과 시각화
plt.figure(figsize=(12, 5))

# 왼쪽 그래프: 액위 변화 (h1, h2)
plt.subplot(1, 2, 1)
plt.plot(t, h1, label='h_1')
plt.plot(t, h2, '--', label='h_2')
plt.xlabel('t(min)')
plt.ylabel('h_1, h_2(m)')
plt.legend(loc='best')
plt.grid(True)

# 오른쪽 그래프: 온도 변화 (T1, T2)
plt.subplot(1, 2, 2)
plt.plot(t, T1, label='T_1')
plt.plot(t, T2, '--', label='T_2')
plt.xlabel('t(min)')
plt.ylabel('T_1, T_2 (deg.C)')
plt.legend(loc='best')
plt.grid(True)

plt.tight_layout()
plt.show()