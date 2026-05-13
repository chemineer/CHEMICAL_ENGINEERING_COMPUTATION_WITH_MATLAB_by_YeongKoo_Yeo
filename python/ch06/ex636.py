import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 파라미터 설정
mumax = 0.33
Cps = 93
Ks = 1.7
kd = 0.01
m = 0.03
Ysc = 1 / 0.08
Ypc = 5.6

# 2. 미분 방정식 정의 (MATLAB의 gf 함수)
def gf(t, C, mumax, Cps, Ks, kd, m, Ysc, Ypc):
    # C[0]=Cc(세포), C[1]=Cs(기질), C[2]=Cp(생성물)
    Cc, Cs, Cp = C
    
    # 속도 식 계산
    rd = kd * Cc  # 사멸 속도[cite: 15]
    rsm = m * Cc  # 유지 대사 속도[cite: 15]
    # 성장 속도 (생성물 저해 효과 포함)[cite: 15]
    rg = mumax * (1 - Cp / Cps)**0.52 * Cc * Cs / (Ks + Cs)
    
    # 농도 변화율[cite: 15]
    dCc = rg - rd
    dCs = -rg * Ysc - rsm
    dCp = Ypc * rg
    
    return [dCc, dCs, dCp]

# 3. 초기 조건 및 시간 구간 설정[cite: 15]
C0 = [1, 250, 0]  # 초기 농도[cite: 15]
tspan = [0, 12]   # 시간 구간 (hr)[cite: 15]

# 4. ODE 풀이 (ode45에 해당)[cite: 15]
sol = solve_ivp(
    gf, 
    tspan, 
    C0, 
    args=(mumax, Cps, Ks, kd, m, Ysc, Ypc), 
    method='RK45', 
    dense_output=True
)

t = sol.t
Cc = sol.y[0, :]
Cs = sol.y[1, :]
Cp = sol.y[2, :]

# 5. 반응 속도 추가 계산 (그래프용)[cite: 15]
rd = kd * Cc
rsm = m * Cc
rg = mumax * (1 - Cp / Cps)**0.52 * Cc * Cs / (Ks + Cs)

# 6. 결과 시각화[cite: 15]
plt.figure(figsize=(12, 10))

# Subplot 1: 세포 농도[cite: 15]
plt.subplot(2, 2, 1)
plt.plot(t, Cc)
plt.xlabel('t(hr)')
plt.ylabel('C_c(g/dm^3)')
plt.grid(True)

# Subplot 2: 기질 및 생성물 농도[cite: 15]
plt.subplot(2, 2, 2)
plt.plot(t, Cs, label='C_s')
plt.plot(t, Cp, '--', label='C_p')
plt.xlabel('t(hr)')
plt.ylabel('C(g/dm^3)')
plt.legend(loc='best')
plt.grid(True)

# Subplot 3: 반응 속도 프로파일[cite: 15]
plt.subplot(2, 2, 3)
plt.plot(t, rg, label='r_g')
plt.plot(t, rsm, '.-', label='r_sm')
plt.plot(t, rd, '--', label='r_d')
plt.xlabel('t(hr)')
plt.ylabel('rates(g/dm^3/hr)')
plt.legend(loc='best')
plt.grid(True)

plt.tight_layout()
plt.show()

# 7. 최종 결과 출력[cite: 15]
print(f'Final concentrations: Ccf = {Cc[-1]:g}, Csf = {Cs[-1]:g}, Cpf = {Cp[-1]:g}')
print(f'Final reaction rates: rgf = {rg[-1]:g}, rsmf = {rsm[-1]:g}, rdf = {rd[-1]:g}')