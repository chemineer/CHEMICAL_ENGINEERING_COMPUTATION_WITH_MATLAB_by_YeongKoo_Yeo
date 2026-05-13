import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 파라미터 설정
Dab = 1.991e-5
T = 328.5
R = 8.314
P = 99.4
Pa0 = 68.4
Tf = 295
z1 = 0
z2 = 0.238
zv = [z1, z2]

# 2. 기초 계산
C = P / (R * T)
x0 = Pa0 / P
critN = 1e-11
errN = 1
Na1 = 3.5e-6
Na2 = 3.6e-6

# 3. 미분 방정식 정의 (dx/dz)
def dxz(z, x, Na, Dab, C):
    # x_A의 거리에 따른 변화율: dx/dz = -(1-x)*Na / (Dab*C)
    return -(1 - x) * Na / (Dab * C)

# 4. 이분법(Bisection Method)을 이용한 Na 추정
Nam = (Na1 + Na2) / 2
xm_final = 0

while errN > critN:
    Nam = (Na1 + Na2) / 2
    
    # 각 Na 후보값에 대해 ODE 풀이
    sol1 = solve_ivp(dxz, zv, [x0], args=(Na1, Dab, C), method='RK45')
    sol2 = solve_ivp(dxz, zv, [x0], args=(Na2, Dab, C), method='RK45')
    solm = solve_ivp(dxz, zv, [x0], args=(Nam, Dab, C), method='RK45')
    
    # 끝점(z2)에서의 xA 값 확인
    x1_end = sol1.y[0, -1]
    x2_end = sol2.y[0, -1]
    xm_end = solm.y[0, -1]
    
    # 이분법 논리 적용
    if x1_end * xm_end < 0:
        Na2 = Nam
    else:
        Na1 = Nam
        
    errN = abs(Na1 - Na2)
    xm_final = xm_end
    z_plot = solm.t
    x_plot = solm.y[0]

# 5. 분석해(Analytic Solution) 계산
xblm = x0 / (np.log(1 / (1 - x0)))
Nanal = Dab * C * x0 / ((z2 - z1) * xblm)

# 6. 결과 출력
print(f'Estimated Nab = {Nam:e}, xA = {xm_final:8.6f}')
print(f'Analytic Nab = {Nanal:e}')

# 7. 시각화
plt.figure(figsize=(8, 5))
plt.plot(z_plot, x_plot)
plt.xlabel('z(m)')
plt.ylabel('x_A')
plt.grid(True)
plt.axis('tight')
plt.title('Concentration Profile of Component A')
plt.show()