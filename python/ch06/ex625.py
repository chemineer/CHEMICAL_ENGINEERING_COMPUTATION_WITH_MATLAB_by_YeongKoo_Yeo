import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from pbrmult import pbrmult  # pbrmult.py 임포트

# 1. 데이터 설정
ka = 100
kc = 1500
Ct0 = 0.2
Ft0 = 20
alpha = 0.0019

# 2. 초기 조건 및 구간 설정
x0 = [10, 10, 0, 0, 1]  # [Fa0, Fb0, Fc0, Fd0, y0]
Wv = [0, 1000]          # W 구간

# 3. ODE 풀이
sol = solve_ivp(
    pbrmult, 
    Wv, 
    x0, 
    args=(ka, kc, Ct0, Ft0, alpha), 
    method='RK45',
    dense_output=True
)

W = sol.t
Fa, Fb, Fc, Fd, y = sol.y

# 4. 선택도 (Selectivity Scd) 계산[cite: 2]
n = len(W)
Scd = np.zeros(n)
for i in range(n):
    if Fd[i] <= 1e-4:
        Scd[i] = 0
    else:
        Scd[i] = Fc[i] / Fd[i]

# 5. 결과 시각화[cite: 2]
plt.figure(figsize=(12, 5))

# Subplot 1: 유량 및 y 변화
plt.subplot(1, 2, 1)
plt.plot(W, Fa, label='F_A', linestyle='-')
plt.plot(W, Fb, label='F_B', linestyle=':')
plt.plot(W, Fc, label='F_C', linestyle='-.')
plt.plot(W, Fd, label='F_D', linestyle='--')
plt.plot(W, y, label='y', linestyle='none', marker='.', markersize=2)
plt.xlabel('W')
plt.ylabel('F_i')
plt.legend()

# Subplot 2: 선택도 Scd
plt.subplot(1, 2, 2)
plt.plot(W, Scd)
plt.xlabel('W')
plt.ylabel('S_{C/D}')

plt.tight_layout()
plt.show()

# 6. 결과 출력[cite: 2]
print("Final molar flow rate of each species:")
print(f" Faf={Fa[-1]:g}, Fbf={Fb[-1]:g}, Fcf={Fc[-1]:g}, Fdf={Fd[-1]:g}")
print(f"Final value of y: yf = {y[-1]:g}")
print(f"Selectivity: Scdf = {Scd[-1]:g}")