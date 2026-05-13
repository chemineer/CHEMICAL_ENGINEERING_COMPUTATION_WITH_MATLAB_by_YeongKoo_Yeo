import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def bf(x1, z, Pt, A, B, C, c, K):
    # z[0] = T, z[1] = L
    T_curr = z[0]
    L_curr = z[1]
    
    x2 = 1 - x1
    x = np.array([x1, x2])
    
    # 증기압 계산 (Antoine 식)
    Pj = 10**(A - B / (T_curr + C))
    
    # 활동도 계수 (Activity coefficients) 계산
    gam = np.zeros(2)
    gam[0] = 10**((1 - x1)**2 * (c[0] + 2 * x1 * (c[1] - c[0])))
    gam[1] = 10**((1 - x2)**2 * (c[1] + 2 * x2 * (c[0] - c[1])))
    
    # 평형 상수 k 계산
    k = gam * Pj / Pt
    
    # 미분 방정식 정의
    dT = K * (1 - k[0] * x1 - k[1] * x2)
    dL = L_curr / (x1 * (k[0] - 1))
    
    return [dT, dL]

# 상수 및 초기값 설정
A = np.array([7.96681, 8.04494])
B = np.array([1668.21, 1554.3])
C = np.array([228, 222.65])
c = np.array([0.3781, 0.6848])
Pt = 760
K = 5e5
x1_span = (0.4, 0.8) # x1의 범위
z0 = [79, 100]       # 초기 온도 및 액체량

# ODE 풀기 (solve_ivp는 기본적으로 RK45 알고리즘 사용)
sol = solve_ivp(
    fun=bf, 
    t_span=x1_span, 
    y0=z0, 
    args=(Pt, A, B, C, c, K),
    method='RK45',
    dense_output=True
)

# 결과 추출
x1_plot = sol.t
T_final = sol.y[0, -1]
L_final = sol.y[1, -1]

print(f"Tfinal = {T_final:g}, Lfinal = {L_final:g}")

# 그래프 출력
plt.figure(figsize=(12, 5))

# 온도 그래프
plt.subplot(1, 2, 1)
plt.plot(sol.t, sol.y[0])
plt.xlabel('x1')
plt.ylabel('T(C)')
plt.grid(True)

# 액체량 그래프
plt.subplot(1, 2, 2)
plt.plot(sol.t, sol.y[1])
plt.xlabel('x1')
plt.ylabel('L(kgmole)')
plt.grid(True)

plt.tight_layout()
plt.show()