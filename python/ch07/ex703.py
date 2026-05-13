import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 파라미터 설정
Dab = 7.39e-6
T = 298.15
P = 1.01325e5
p0 = 133.32
r0 = 0.003
xspan = [0, 20]
critN = 1e-16
errN = 1
Na1 = 1e-12
Na2 = 2e-12  # 초기 추측치

# 2. 미분 방정식 정의 (spf 함수)
def spf(x, z, T, P, Dab, r0):
    # z[0] = Na * r^2 (상수), z[1] = 부분압 p
    R = 8314.34
    # dz/dx 계산
    # dz[0]/dx = 0 (연속 방정식에 의해 Na*r^2은 일정함)
    # dz[1]/dx = dp/dx 식 적용[cite: 13]
    dzdx = [0, -R * T * z[0] * (1 - z[1] / P) / (Dab * (x + r0)**2)]
    return dzdx

# 3. 이분법을 이용한 Na 추정[cite: 13]
Nam = (Na1 + Na2) / 2
sol_m = None

while errN > critN:
    Nam = (Na1 + Na2) / 2
    
    # 각 Na 후보값(z[0]에 해당)에 대해 ODE 풀이[cite: 13]
    # xspan[0]에서 p=p0인 초기 조건 적용[cite: 13]
    sol1 = solve_ivp(spf, xspan, [Na1, p0], args=(T, P, Dab, r0), method='RK45')
    sol2 = solve_ivp(spf, xspan, [Na2, p0], args=(T, P, Dab, r0), method='RK45')
    sol_m = solve_ivp(spf, xspan, [Nam, p0], args=(T, P, Dab, r0), method='RK45')
    
    # 무한대 거리(끝점)에서의 부분압 부호 확인을 통한 범위 수정[cite: 13]
    p_end1 = sol1.y[1, -1]
    p_endm = sol_m.y[1, -1]
    
    if p_end1 * p_endm < 0:
        Na2 = Nam
    else:
        Na1 = Nam
        
    errN = abs(Na1 - Na2)

# 4. 결과 데이터 계산 및 출력[cite: 13]
x = sol_m.t
zm = sol_m.y
r = r0 + x
# Na = (Na * r^2) / r^2[cite: 13]
Na_profile = zm[0, :] / r**2
pf = zm[1, -1]

print(f'Flux at r=r0: {Na_profile[0]:.5e}, partial pressure at r=inf: {pf:.5f}')

# 5. 시각화[cite: 13]
plt.figure(figsize=(8, 5))
# r0 근처의 프로파일 확인 (원본 코드의 r(1:25) 반영)[cite: 13]
plt.plot(r[:25], Na_profile[:25])
plt.xlabel('r(m)')
plt.ylabel('N_A')
plt.grid(True)
plt.title('Molar Flux Profile (Spherical Diffusion)')
plt.show()