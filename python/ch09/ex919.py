import matplotlib.pyplot as plt
import numpy as np
from stepnp import stepnp  # stepnp.py에서 함수 임포트

# 1. 기본 데이터 설정
tau = 2
zeta = 0.25
Kp = 3
Kc = 5
h = ['-', ':', '-.', '--'] # 선 스타일
t0 = 0
delt = 0.1
fint = 20
ms = 1

plt.figure(figsize=(12, 5))

# --- 좌측 그래프: Constant reset time (tauI = 1) ---
tauI = 1
tauD_list = [0.5, 1, 5, 10]

plt.subplot(1, 2, 1)
for i, tauD in enumerate(tauD_list):
    # 분자(num) 및 분모(den) 계수 계산
    num = [Kc * Kp * tauI * tauD, Kc * Kp * tauI, Kc * Kp]
    d1 = tauI * (tau**2)
    d2 = 2 * tauI * tau * zeta + Kc * Kp * tauI * tauD
    d3 = tauI * (1 + Kc * Kp)
    d4 = Kc * Kp
    den = [d1, d2, d3, d4]
    
    # stepnp 호출
    y, t = stepnp(num, den, t0, delt, fint, ms)
    plt.plot(t, y, linestyle=h[i], label=f'$\\tau_D$={tauD}')

# 설정 및 데코레이션
st = np.ones_like(t) # 목표값(1) 표시
plt.plot(t, st, 'k-', linewidth=0.5)
plt.xlabel('Time(min)')
plt.ylabel('Output, y(t)')
plt.title('PID control(Kc=5, \\tau_I=1)')
plt.legend(loc='best')
plt.grid(True)

# --- 우측 그래프: Constant derivative time (tauD = 0.5) ---
tauD = 0.5
tauI_list = [0.5, 1, 5, 10]

plt.subplot(1, 2, 2)
for i, tauI in enumerate(tauI_list):
    # 분자(num) 및 분모(den) 계수 계산
    num = [Kc * Kp * tauI * tauD, Kc * Kp * tauI, Kc * Kp]
    d1 = tauI * (tau**2)
    d2 = 2 * tauI * tau * zeta + Kc * Kp * tauI * tauD
    d3 = tauI * (1 + Kc * Kp)
    d4 = Kc * Kp
    den = [d1, d2, d3, d4]
    
    # stepnp 호출
    y, t = stepnp(num, den, t0, delt, fint, ms)
    plt.plot(t, y, linestyle=h[i], label=f'$\\tau_I$={tauI}')

# 설정 및 데코레이션
plt.plot(t, st, 'k-', linewidth=0.5)
plt.xlabel('Time(min)')
plt.ylabel('Output, y(t)')
plt.title('PID control(Kc=5, \\tau_D=0.5)')
plt.legend(loc='best')
plt.grid(True)

plt.tight_layout()
plt.show()