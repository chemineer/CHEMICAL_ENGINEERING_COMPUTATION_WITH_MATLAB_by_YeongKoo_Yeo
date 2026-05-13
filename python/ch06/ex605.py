import numpy as np
from scipy.optimize import fsolve

# --- 데이터 및 상수 설정 ---
P = 2           # 압력 (atm)
T = 340         # 온도 (K)
R = 0.082       # 기체 상수 (L·atm/mol·K)
Kc = 0.1        # 평형 상수[cite: 9]
ya0 = 1         # 초기 몰 분율[cite: 9]
epsilon = 1     # 팽창 계수 (Expansion factor)[cite: 9]

# 초기 농도 계산: Ca0 = ya0 * P / (R * T)[cite: 9]
Ca0 = ya0 * P / (R * T)

# --- 평형 방정식 정의 ---[cite: 9]
# fF = @(x) 4*Ca0*x^2 - Kc*(1-x)*(1+epsilon*x)[cite: 9]
def fF(x):
    return 4 * Ca0 * x**2 - Kc * (1 - x) * (1 + epsilon * x)

# --- 수치 해법 실행 ---
# 초기 추정값 x0 = 0.5[cite: 9]
x0 = 0.5

# fsolve를 사용하여 평형 전환율 Xe를 구합니다.
# MATLAB의 fzero와 대응되는 라이브러리입니다.[cite: 9]
solution = fsolve(fF, x0)
Xe = solution[0]

# --- 결과 출력 ---
print(f"Initial Concentration (Ca0): {Ca0:.6f}")
print(f"Equilibrium Conversion (Xe): {Xe:.6f}")