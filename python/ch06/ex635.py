import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 파라미터 설정
k1 = 2
k2 = 3
Ct0 = 0.8
Fb0 = 4
Vt = 50
Fai = 4
Rb = Fb0 / Vt # 성분 B의 공급 속도

# 2. 미분 방정식 정의 (MATLAB의 mbf 함수)
def mbf(V, x, k1, k2, Ct0, Rb):
    # x[0]=Fa, x[1]=Fb, x[2]=Fd, x[3]=Fu
    # 농도 계산: Ci = Ct0 * (Fi / sum(Fi))[cite: 14]
    sum_x = np.sum(x)
    if sum_x == 0:
        C = np.zeros(4)
    else:
        C = Ct0 * x / sum_x
    
    # 반응 속도 및 유량 변화율 정의[cite: 14]
    # Fa 변화: -k1*Ca^2*Cb - k2*Ca*Cb^2
    dfa = -k1 * C[0]**2 * C[1] - k2 * C[0] * C[1]**2
    # Fb 변화: Fa 변화량 + 공급량(Rb)
    dfb = -k1 * C[0]**2 * C[1] - k2 * C[0] * C[1]**2 + Rb
    # Fd 변화 (원하는 생성물)
    dfd = k1 * C[0]**2 * C[1]
    # Fu 변화 (원하지 않는 생성물)
    dfu = k2 * C[0] * C[1]**2
    
    return [dfa, dfb, dfd, dfu]

# 3. 초기 조건 및 ODE 풀이[cite: 14]
x0 = [Fai, 0, 0, 0] # 초기 몰 유량[cite: 14]
Vspan = [0, Vt]     # 부피 구간[cite: 14]

sol = solve_ivp(mbf, Vspan, x0, args=(k1, k2, Ct0, Rb), method='RK45', dense_output=True)

V = sol.t
Fa, Fb, Fd, Fu = sol.y

# 4. 선택도 (Sdu) 계산[cite: 14]
n = len(V)
Sdu = np.zeros(n)
for i in range(n):
    # 분모가 0이 되는 것을 방지하기 위한 조건문[cite: 14]
    if Fu[i] <= 1e-6:
        Sdu[i] = 0
    else:
        Sdu[i] = Fd[i] / Fu[i]

# 5. 결과 시각화[cite: 14]
plt.figure(figsize=(12, 5))

# Subplot 1: 유량 프로파일[cite: 14]
plt.subplot(1, 2, 1)
plt.plot(V, Fa, label='F$_A$', linestyle='-')
plt.plot(V, Fb, label='F$_B$', linestyle=':')
plt.plot(V, Fd, label='F$_D$', linestyle='-.', marker='.', markevery=5)
plt.plot(V, Fu, label='F$_U$', linestyle='--')
plt.xlabel('V(dm$^3$)')
plt.ylabel('F$_i$(mol/s)')
plt.legend(loc='best')
plt.grid(True)

# Subplot 2: 선택도 변화[cite: 14]
plt.subplot(1, 2, 2)
plt.plot(V, Sdu)
plt.xlabel('V(dm$^3$)')
plt.ylabel('S$_{D/U}$')
plt.grid(True)

plt.tight_layout()
plt.show()

# 6. 최종 결과 출력[cite: 14]
print('Final molar flow rates:')
print(f'Faf = {Fa[-1]:g},  Fbf = {Fb[-1]:g},  Fdf = {Fd[-1]:g},  Fuf = {Fu[-1]:g}')
print(f'Overall selectivity: Sduf = {Sdu[-1]:g}')