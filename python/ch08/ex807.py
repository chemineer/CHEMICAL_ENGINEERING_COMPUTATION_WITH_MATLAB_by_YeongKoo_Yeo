import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 초기 데이터 설정
T0 = 260; Ta = 60; R1 = 0.0833; R2 = 0.25; B = 0.5/12; k = 133; h = 3
qa = 200; qb = 800  # 열유속 추측 범위
feps = 1e-3; ferr = 1e3; iter_count = 1

# 2. 미분 방정식 정의: dy/dr
def funQL(r, y, Ta, k, h, B):
    T_val = y[0]
    y2_val = y[1]
    
    q_val = y2_val / r
    dTdr = -q_val / k
    dy2dr = -h * r * (T_val - Ta) / B
    
    return [dTdr, dy2dr]

# 3. 이분법 루프 (경계 조건 만족 찾기)
r_eval = np.linspace(R1, R2, 100)
ym_final = None

while ferr > feps:
    qm = (qa + qb) / 2
    
    # 각 추측값에 대해 ODE 풀이
    # sol.y[0] = T, sol.y[1] = r * q
    sol_a = solve_ivp(funQL, [R1, R2], [T0, qa], args=(Ta, k, h, B), t_eval=[R2])
    sol_m = solve_ivp(funQL, [R1, R2], [T0, qm], args=(Ta, k, h, B), t_eval=r_eval)
    
    # r=R2에서의 함수 값 계산 (Boundary condition: q - h(T-Ta) = 0 형태)
    # MATLAB의 y(end, 2)는 파이썬에서 sol.y[1, -1]에 대응
    fa = sol_a.y[1, -1] - R2 * h * (sol_a.y[0, -1] - Ta)
    fm = sol_m.y[1, -1] - R2 * h * (sol_m.y[0, -1] - Ta)
    
    if (fa * fm) < 0:
        qb = qm
    else:
        qa = qm
        
    ferr = abs(qb - qa)
    iter_count += 1
    ym_final = sol_m

# 4. 결과 데이터 가공
r_plot = ym_final.t
T_plot = ym_final.y[0]
y2_plot = ym_final.y[1]
qrate_plot = 2 * np.pi * y2_plot * B # 열전달률 (Btu/hr)

# 5. 시각화
plt.figure(figsize=(12, 5))

# 왼쪽 그래프: 열전달률 (qrate)
plt.subplot(1, 2, 1)
plt.plot(r_plot, qrate_plot)
plt.xlabel('r(ft)')
plt.ylabel('q(Btu/hr)')
plt.grid(True)
plt.autoscale(enable=True, axis='x', tight=True)

# 오른쪽 그래프: 온도 (T)
plt.subplot(1, 2, 2)
plt.plot(r_plot, T_plot)
plt.xlabel('r(ft)')
plt.ylabel('T(F)')
plt.grid(True)

plt.tight_layout()
plt.show()

# 6. 결과 출력
print(f"Number of iterations: {iter_count}")
print(f"heat transfer rate at r=R2 (Btu/h): {qrate_plot[-1]:.6g}")
print(f"T at r=R2 (F): {T_plot[-1]:.6g}")