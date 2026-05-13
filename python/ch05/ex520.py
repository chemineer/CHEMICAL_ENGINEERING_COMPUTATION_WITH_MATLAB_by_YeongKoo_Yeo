import numpy as np
from scipy.optimize import fsolve

# 데이터 설정 (두 코드 동일)
G = np.float64(3900 / 3600)  # 질량 유량 (kg/s)
d1 = 0.10226     # 상류 배관 직경 (m)
d2 = 0.035       # 오리피스 직경 (m)
p1 = np.float64(1.2e6)       # 상류 압력 (Pa)
Fa = 1           # 열팽창 계수
rho = 10.25      # 밀도 (kg/m^3)
mu = np.float64(1.3e-5)      # 점도 (Pa*s)
gam = 1.3        # 비열비

# 계산
beta = d2 / d1
Ao = np.pi * d2**2 / 4  # 오리피스 면적
Cd = 0.6274 - 0.2354 * beta + 0.7858 * beta**2

# 비선형 방정식 정의
def objective(dp):
    # dp가 0 이하일 경우를 대비하여 방어 코드 추가
    if dp <= 0: return 1e6
    term1 = 1 - (0.41 + 0.35 * beta**4) * dp / (gam * p1)
    term2 = (np.sqrt((rho * (1 - beta**4)) / (2 * dp))) * (G / (rho * Ao * Cd * Fa))
    return term1 - term2

# 수정된 부분: MATLAB과 동일한 초기값 사용
x0 = 116 
dp_solution, info, ier, mesg = fsolve(objective, x0, full_output=True, xtol=1e-12)
dp = dp_solution[0]

# 팽창 계수(Y) 계산
Y = 1 - (0.41 + 0.35 * beta**4) * dp / (gam * p1)

# 결과 출력
print(f"Expansion factor = {Y:.6f}")
print(f"Pressure drop = {dp / 1000:.6f} kPa")
print(ier)