import numpy as np
import matplotlib.pyplot as plt

# 1. 초기 설정 및 데이터 정의
T = 200                # 온도 (K)
R = 0.08206            # 기체 상수 (L·atm/mol·K)

# 압력 데이터 P(atm) 및 밀도 데이터 rho(mol/l)
P = np.array([3.2, 6.0, 9.0, 12.0, 14.0, 17.0, 19.0, 21.0])
rho = np.array([0.1995, 0.3825, 0.5895, 0.8108, 0.9685, 1.2257, 1.4159, 1.6281])

# 2. 압축 인자 Z 계산 및 다항식 적합 (Polynomial Fitting)
# Z = P / (R * rho * T)
Z = P / (R * rho * T)

# rho와 Z의 관계를 3차 다항식으로 근사
Vc = np.polyfit(rho, Z, 3)
print(f"Polynomial Coefficients (Z vs rho): \n{Vc}")

# 3. 피팅 결과 계산
# x축 범위 설정 (0.15부터 1.8까지 0.01 간격)
x = np.arange(0.15, 1.81, 0.01)

# 근사된 다항식을 이용해 Z값(Zcal) 계산
Zcal = np.polyval(Vc, x)

# P = R * Z * rho * T 관계식을 이용해 계산된 압력(Pcal) 도출
Pcal = R * Zcal * x * T

# 4. 시각화
plt.figure(figsize=(10, 6))
plt.plot(x, Pcal, label='Fitting curve') # 피팅 곡선
plt.plot(rho, P, 'o', label='Data')       # 원본 데이터 점
plt.grid(True)
plt.xlabel(r'$\rho$(mol/liter)')          # LaTeX 스타일 수식 적용
plt.ylabel('P(atm)')
plt.legend(loc='best')
plt.title('N2 Density Fitting (Pressure vs Density)')
plt.show()