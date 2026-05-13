import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 미분 방정식 정의 (MATLAB의 pbconv 함수)
def pbconv(w, z, kp, Fa0, epsilon, alpha):
    # x = z[0], y = z[1]
    dx_dw = (kp * (1 - z[0]) * z[1]) / (Fa0 * (1 + epsilon * z[0]))
    dy_dw = -alpha * (1 + epsilon * z[0]) / (2 * z[1])
    return [dx_dw, dy_dw]

# 2. 데이터 및 초기 파라미터 설정
k = 0.00392
Fa0 = 0.1362
Fb0 = 0.068
P0 = 10
alpha = 0.0367

Fi = Fb0 * (79 / 21)
Ft0 = Fa0 + Fb0 + Fi
ya0 = Fa0 / Ft0
Pa0 = ya0 * P0
kp = k * Pa0 * (0.5)**(2/3)
delta = -0.5
epsilon = ya0 * delta

Wfa = 15.0  # 초기 추측치 A
Wfb = 25.0  # 초기 추측치 B
z0 = [0, 1] # 초기 조건 [X=0, y=1][cite: 3]
Xf = 0.6    # 목표 전환율[cite: 3]

# 3. 이분법(Bisection Method)을 이용한 촉매 무게 결정[cite: 3]
while abs(Wfa - Wfb) >= 1e-3:
    Wfm = (Wfa + Wfb) / 2
    
    # 각 구간별 ODE 풀이
    sol_a = solve_ivp(pbconv, [0, Wfa], z0, args=(kp, Fa0, epsilon, alpha))
    sol_m = solve_ivp(pbconv, [0, Wfm], z0, args=(kp, Fa0, epsilon, alpha))
    
    Xa = sol_a.y[0, -1] # Wfa일 때의 최종 전환율
    Xm = sol_m.y[0, -1] # Wfm일 때의 최종 전환율
    
    # 이분법 조건 체크[cite: 3]
    if (Xa - Xf) * (Xm - Xf) < 0:
        Wfb = Wfm
    else:
        Wfa = Wfm

# 최종 결정된 구간으로 다시 계산
sol_final = solve_ivp(pbconv, [0, Wfm], z0, args=(kp, Fa0, epsilon, alpha), dense_output=True)
W_plot = np.linspace(0, Wfm, 100)
z_plot = sol_final.sol(W_plot)

X = z_plot[0]
y = z_plot[1]
Wm = W_plot

# 4. 결과 출력 및 추가 계산[cite: 3]
print(f'Catalyst weight = {Wfm:g}, Conversion = {X[-1]:g}')

fm = (1 + epsilon * X) / y  # 부피 유량 비율[cite: 3]
rp = -kp * (1 - X) * y / (1 + epsilon * X)  # 반응 속도[cite: 3]

# 5. 시각화[cite: 3]
plt.figure(figsize=(12, 5))

# Subplot 1: X, y, f 변화
plt.subplot(1, 2, 1)
plt.plot(Wm, X, label='X', linestyle='-')
plt.plot(Wm, y, label='y', linestyle=':')
plt.plot(Wm, fm, label='f', linestyle='--')
plt.xlabel('W(kg)')
plt.ylabel('Values')
plt.legend()
plt.grid(True)

# Subplot 2: 반응 속도 -rA
plt.subplot(1, 2, 2)
plt.plot(Wm, -rp)
plt.xlabel('W(kg)')
plt.ylabel('-r_A')
plt.grid(True)

plt.tight_layout()
plt.show()