import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from catLH import catLH

# 2. 데이터 및 파라미터 설정[cite: 21]
n = 50                 # 격자 분할 수
R = 0.2                # 촉매 입자 반지름 (cm)
k = 100                # 반응 속도 상수
De = 0.25              # 유효 확산 계수 (cm^2/s)
Kr = 1e9               # 흡착 평형 상수 (cm^3/mol)
h = R / n              # 격자 간격
Ca0_list = 5e-5 * np.array([1, 2, 3, 4]) # 표면 농도 리스트[cite: 21]

# 3. 결과 저장을 위한 설정
x = np.linspace(0, R, n + 1) / R # 무차원 거리 (r/R)[cite: 21]
results = []

# 4. 각 표면 농도에 대하여 비선형 방정식 풀이[cite: 21]
for Cas in Ca0_list:
    # 초기 추정값: 표면 농도로 균일하다고 가정[cite: 21]
    y0 = Cas * np.ones(n + 1)
    
    # fsolve를 이용한 수치 해 도출
    y_sol = fsolve(catLH, y0, args=(De, n, h, k, Kr, Cas))
    
    # 무차원 농도 (Ca/Cas) 저장[cite: 21]
    results.append(y_sol / Cas)

# 5. 시각화[cite: 21]
plt.figure(figsize=(8, 6))
line_styles = ['-', ':', '--', '-.']
labels = [f'$C_{{A0}} = {val*1e5:.0f} \\times 10^{{-5}}$ mol/cm$^3$' for val in Ca0_list]

for i in range(len(Ca0_list)):
    plt.plot(x, results[i], line_styles[i], label=labels[i])

plt.xlabel('r/R (Dimensionless Radius)')
plt.ylabel('$C_A/C_{A0}$ (Dimensionless Concentration)')
plt.title('Diffusion and Reaction in Catalyst Pellet (L-H Kinetics)')
plt.legend(loc='best')
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()