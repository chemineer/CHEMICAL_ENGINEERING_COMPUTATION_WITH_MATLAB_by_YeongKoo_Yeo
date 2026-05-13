import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 파라미터 및 초기 설정
Cas = 3e-5
Ct = 4e-5
L = 0.2
De = 0.01
k = 8e4
Kc = 6e5
zspan = [0, L]
critN = 1e-10
errN = 1
Na1 = 2e-6
Na2 = 5e-6

# 2. 미분 방정식 정의 (crf 함수)
def crf(z, x, Ct, k, Kc, De):
    # x[0] = Na, x[1] = Ca
    Na = x[0]
    Ca = x[1]
    
    # dx/dz 계산
    # dNa/dz = -k * (Ca^2 - (Ct - Ca)/Kc)
    dNa_dz = -k * (Ca**2 - (Ct - Ca) / Kc)
    # dCa/dz = Na * (Ca/Ct - 2) / (2 * De)
    dCa_dz = Na * (Ca / Ct - 2) / (2 * De)
    
    return [dNa_dz, dCa_dz]

# 3. 이분법을 이용한 초기 플럭스(Na) 추정[cite: 13]
Nam = (Na1 + Na2) / 2
sol_m = None

while errN > critN:
    Nam = (Na1 + Na2) / 2
    
    # 각 Na 후보값에 대해 ODE 풀이[cite: 13]
    # 초기 조건: z=0에서 Na=후보값, Ca=Cas[cite: 13]
    sol1 = solve_ivp(crf, zspan, [Na1, Cas], args=(Ct, k, Kc, De), method='RK45')
    solm = solve_ivp(crf, zspan, [Nam, Cas], args=(Ct, k, Kc, De), method='RK45')
    
    # 끝점(z=L)에서 Na=0 조건을 만족하는지 확인[cite: 13]
    if sol1.y[0, -1] * solm.y[0, -1] < 0:
        Na2 = Nam
    else:
        Na1 = Nam
        
    errN = abs(Na1 - Na2)
    sol_m = solm

# 4. 결과 계산 및 출력[cite: 13]
z = sol_m.t
Na_profile = sol_m.y[0, :]
Ca_profile = sol_m.y[1, :]

# 유효 계수(Effectiveness factor) 계산[cite: 13]
ras = Cas**2 - (Ct - Cas) / Kc
effc = Na_profile[0] / (L * k * ras)

print(f'Effectiveness factor: {effc:7.5f}')
print(f'Molar flux of A at z=0: {Na_profile[0]:12.8f}')

# 5. 시각화 (Subplot 구성)[cite: 13]
plt.figure(figsize=(12, 5))

# 농도 프로파일[cite: 13]
plt.subplot(1, 2, 1)
plt.plot(z, Ca_profile)
plt.xlabel('z(cm)')
plt.ylabel('C_A(gmol/cm^3)')
plt.grid(True)

# 플럭스 프로파일[cite: 13]
plt.subplot(1, 2, 2)
plt.plot(z, Na_profile)
plt.xlabel('z(cm)')
plt.ylabel('N_A(gmol/s/cm^2)')
plt.grid(True)

plt.tight_layout()
plt.show()