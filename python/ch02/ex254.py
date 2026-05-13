import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 파라미터 설정
A = 0.17142
C = 205.74
F = 8000
Hg = 320
Ht = 266.67
Hw = 1.6
Pe = 0.1
Te = 600
Tw = 720

# 2. 미분 방정식 시스템 정의
def flupb_system(t, y):
    # y[0]: P, y[1]: T, y[2]: Pp, y[3]: Tp
    P, T, Pp, Tp = y
    
    # 속도 상수 및 반응 항 계산 (지수 함수 포함)
    reaction_term = 6e-4 * np.exp(20.7 - 15000 / Tp)
    
    # 각 상태 변수에 대한 변화율 (dy/dt)
    dP = Pe - P + Hg * (Pp - P)
    dT = Te - T + Ht * (Tp - T) + Hw * (Tw - T)
    dPp = Hg * (P - Pp * (1 + reaction_term)) / A
    dTp = Ht * ((T - Tp) + F * reaction_term * Pp) / C
    
    return [dP, dT, dPp, dTp]

# 3. 초기 조건 및 시간 설정
y0 = [0.1, 600, 0, 761]
t_span = (0, 1500)
t_eval = np.linspace(0, 1500, 1000) # 부드러운 그래프를 위해 지점 설정

# 4. ODE 풀이 (Stiff 시스템이므로 BDF 방식 사용)
sol = solve_ivp(flupb_system, t_span, y0, method='BDF', t_eval=t_eval)

# 5. 결과 시각화
plt.figure(figsize=(12, 5))

# 왼쪽 그래프: 압력 프로파일 (P, Pp)
plt.subplot(1, 2, 1)
plt.plot(sol.t, sol.y[0], label='P')
plt.plot(sol.t, sol.y[2], ':', label='Pp')
plt.grid(True)
plt.xlabel('τ')
plt.ylabel('Pressure (atm)')
plt.legend()
plt.title('Pressure Profiles')

# 오른쪽 그래프: 온도 프로파일 (T, Tp)
plt.subplot(1, 2, 2)
plt.plot(sol.t, sol.y[1], label='T')
plt.plot(sol.t, sol.y[3], ':', label='Tp')
plt.grid(True)
plt.xlabel('τ')
plt.ylabel('Temperature (R)')
plt.legend()
plt.title('Temperature Profiles')

plt.tight_layout()
plt.show()