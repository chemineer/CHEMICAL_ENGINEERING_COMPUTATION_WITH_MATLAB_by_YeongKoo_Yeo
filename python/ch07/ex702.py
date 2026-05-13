import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 파라미터 설정
Ca1, Cb1, Cc1 = 2.229e-4, 0, 7.208e-3
Ca2, Cb2, Cc2 = 0, 2.701e-3, 4.73e-3
D1, D2, D3 = 1.075e-4, 1.245e-4, 1.47e-4  # Dab, Dbc, Dac
P = 0.2
T = 328
R = 82.057e-3
Ct = P / (R * T)
L = 0.001
Nb = -4.143e-4
Nc = 0

# 2. 초기 추정치 및 반복문 설정
Na_initial = -D3 * (Ca2 - Ca1) / L
zspan = [0, L]
c0 = [Ca1, Cb1, Cc1]
critN = 1e-10
errA = 1
Na1 = Na_initial / 2
Na2 = 2 * Na_initial
iter_count = 1
Nam = (Na1 + Na2) / 2

# 3. 미분 방정식 정의 (Maxwell-Stefan)
def mf(z, c, D1, D2, D3, Na, Nb, Nc, Ct):
    # c[0]=Ca, c[1]=Cb, c[2]=Cc
    xa = c[0] / Ct
    xb = c[1] / Ct
    xc = c[2] / Ct
    
    # dc/dz 계산
    dca_dz = (xa * Nb - xb * Na) / D1 + (xa * Nc - xc * Na) / D3
    dcb_dz = (xb * Na - xa * Nb) / D1 + (xb * Nc - xc * Nb) / D2
    dcc_dz = (xc * Na - xa * Nc) / D3 + (xc * Nb - xb * Nc) / D2
    
    return [dca_dz, dcb_dz, dcc_dz]

# 4. 이분법을 이용한 Na 추정[cite: 13]
while errA > critN:
    Nam = (Na1 + Na2) / 2
    
    # ODE 풀이 (MATLAB의 ode45에 해당하는 RK45 사용)[cite: 13]
    sol1 = solve_ivp(mf, zspan, c0, args=(D1, D2, D3, Na1, Nb, Nc, Ct), method='RK45')
    solm = solve_ivp(mf, zspan, c0, args=(D1, D2, D3, Nam, Nb, Nc, Ct), method='RK45')
    
    # 끝점(z=L)에서의 Ca 농도(Ca2=0) 만족 여부 확인[cite: 13]
    # sol.y[0, -1]은 마지막 지점의 Ca 값
    if sol1.y[0, -1] * solm.y[0, -1] < 0:
        Na2 = Nam
    else:
        Na1 = Nam
        
    errA = abs(Na1 - Na2)
    iter_count += 1

# 5. 결과 데이터 정리 및 시각화[cite: 13]
# 최종 수렴한 Nam으로 다시 풀기 (그래프용)
sol_final = solve_ivp(mf, zspan, c0, args=(D1, D2, D3, Nam, Nb, Nc, Ct), 
                      method='RK45', t_eval=np.linspace(0, L, 100))

z = sol_final.t
c = sol_final.y
xa = c[0, :] / Ct
xb = c[1, :] / Ct
xc = c[2, :] / Ct

# 결과 출력[cite: 13]
print(f"Iterations: {iter_count}")
print(f"Final Na: {Nam:e}")

# 그래프 그리기[cite: 13]
plt.figure(figsize=(8, 6))
plt.plot(z, xa, label='x_A', linestyle='-')
plt.plot(z, xb, label='x_B', linestyle=':')
plt.plot(z, xc, label='x_C', linestyle='-.', marker='.', markevery=10)

plt.xlabel('Distance z(m)')
plt.ylabel('Mole fraction')
plt.legend()
plt.grid(True)
plt.show()