import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 초기 데이터 설정
T1 = 288.15
R1 = 0.004
k = 5
Ir = 400 / (np.pi * R1**2)

# 이분법 초기값 및 설정
T0a = 300
T0b = 400
Teps = 1e-3
Terr = 1e3
iter_count = 1

# 2. 미분 방정식 정의 (dz/dr)
def funQT(r, z, k, R1, Ir):
    T_val = z[0]
    z2_val = z[1]
    
    # r=0에서의 특이점 방지
    if r > 0:
        Qr = z2_val / r
    else:
        Qr = 0
        
    dz1 = -Qr / k
    # 온도에 의존하는 저항/발열 항 계산
    dz2 = (Ir**2 * r) / (1.4e5 * np.exp(0.0035 * T_val))
    
    return [dz1, dz2]

# 3. 이분법 루프 (Boundary Value Problem 해결)
zm_final = None
r_eval = np.linspace(0, R1, 100)

while Terr > Teps:
    T0m = (T0a + T0b) / 2
    
    # 각각의 추측 온도에 대해 ODE 풀이
    # MATLAB의 [T0a 0]은 r=0에서의 [온도, z2값] 초기 조건을 의미함
    sol_a = solve_ivp(funQT, [0, R1], [T0a, 0], args=(k, R1, Ir), t_eval=[R1])
    sol_m = solve_ivp(funQT, [0, R1], [T0m, 0], args=(k, R1, Ir), t_eval=r_eval)
    
    # r=R1(마지막 지점)에서의 온도 값 추출
    ta_end = sol_a.y[0, -1]
    tm_end = sol_m.y[0, -1]
    
    # 이분법 조건 체크 (목표값 T1과의 차이)
    if (tm_end - T1) > 0:
        T0b = T0m
    else:
        T0a = T0m
        
    Terr = abs(T0b - T0a)
    iter_count += 1
    zm_final = sol_m # 최종 결과 저장을 위해 업데이트

# 4. 결과 데이터 정리
r_plot = zm_final.t
T_plot = zm_final.y[0]
# Qr = z(2)/r 계산 (r=0 제외 처리)
Qr_plot = np.zeros_like(r_plot)
Qr_plot[1:] = zm_final.y[1, 1:] / r_plot[1:]

# 5. 시각화
plt.figure(figsize=(12, 5))

# 왼쪽 그래프: 열유속 (Qr)
plt.subplot(1, 2, 1)
plt.plot(r_plot, Qr_plot)
plt.xlabel('r(m)')
plt.ylabel('Qr(W/m^2)')
plt.grid(True)

# 오른쪽 그래프: 온도 (T)
plt.subplot(1, 2, 2)
plt.plot(r_plot, T_plot)
plt.xlabel('r(m)')
plt.ylabel('T(K)')
plt.grid(True)

plt.tight_layout()
plt.show()

# 6. 결과 출력
print(f"Number of iterations: {iter_count}")
print(f"Qr at r=R1 (W/m^2): {Qr_plot[-1]:.6g}")
print(f"T at r=0 (K): {T_plot[0]:.6g}")