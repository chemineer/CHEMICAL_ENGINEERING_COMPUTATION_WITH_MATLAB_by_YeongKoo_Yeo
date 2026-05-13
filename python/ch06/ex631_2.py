import numpy as np
from scipy.integrate import solve_ivp
from adfun import adfun  # adfun.py 모듈 임포트
import matplotlib.pyplot as plt

#  convtemp.m 
# 1. 데이터 및 초기 설정
P = 162
Vspan = [0, 4]
# 질소 유량 배열
FN2_list = np.array([28.3, 18.3, 8.3, 3.3, 0.0])
nF = len(FN2_list)

# 결과를 저장할 빈 배열
xc = np.zeros(nF)
Tr = np.zeros(nF)

# 2. 루프를 통한 조건별 계산[cite: 9]
for i in range(nF):
    fn2 = FN2_list[i]
    # 초기 조건: [FA0, FB0, FC0, T0][cite: 9]
    # FA0는 전체 유량(38.3)에서 질소 유량을 뺀 값[cite: 9]
    X0 = [38.3 - fn2, 0, 0, 1150]
    pf = [P, fn2]
    
    # ODE 풀이[cite: 9]
    sol = solve_ivp(adfun, Vspan, X0, args=(pf,), method='RK45')
    
    # 최종 결과 저장[cite: 9]
    # 전환율 xc = (FA0 - FA_final) / FA0[cite: 9]
    xc[i] = (X0[0] - sol.y[0, -1]) / X0[0]
    # 최종 온도 Tr[cite: 9]
    Tr[i] = sol.y[3, -1]

# 3. 결과 출력[cite: 9]
print("Final Conversions (xc):", xc)
print("Final Temperatures (Tr):", Tr)


#  finxtemp.m 
# 1. 특정 데이터 설정[cite: 9]
P = 162
Vspan = [0, 4]
FN2 = 28.3
X0 = [38.3 - FN2, 0, 0, 1150]
pf = [P, FN2]

# 2. ODE 풀이[cite: 9]
# 그래프를 위해 dense_output=True 사용[cite: 9]
sol = solve_ivp(adfun, Vspan, X0, args=(pf,), method='RK45', dense_output=True)

V = np.linspace(0, 4, 100)
X_results = sol.sol(V)

# 3. 데이터 가공[cite: 9]
# 전환율 xc 계산: (FA0 - FA(V)) / FA0[cite: 9]
xc_profile = (X0[0] - X_results[0, :]) / X0[0]
temp_profile = X_results[3, :] # 온도 프로파일[cite: 9]

# 4. 시각화[cite: 9]
plt.figure(figsize=(12, 5))

# Subplot 1: 온도 프로파일[cite: 9]
plt.subplot(1, 2, 1)
plt.plot(V, temp_profile)
plt.xlabel('Reactor volume(m$^3$)')
plt.ylabel('Temperature(K)')
plt.grid(True)

# Subplot 2: 전환율 프로파일[cite: 9]
plt.subplot(1, 2, 2)
plt.plot(V, xc_profile)
plt.xlabel('Reactor volume(m$^3$)')
plt.ylabel('Conversion')
plt.grid(True)

plt.tight_layout()
plt.show()