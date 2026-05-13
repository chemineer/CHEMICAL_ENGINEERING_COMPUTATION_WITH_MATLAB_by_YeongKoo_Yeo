import numpy as np
from scipy.optimize import fsolve

# 1. 입력 데이터 및 상수 설정
R = 0.082054       # 기체 상수 (L·atm/mol·K)
Tc = 419.6         # 임계 온도 (Critical Temperature, K)
Pc = 396.743       # 임계 압력 (Critical Pressure, atm)
T = 415            # 현재 온도 (K)
P = 207.2538       # 현재 압력 (atm)

# 2. SRK 매개변수 a, b 계산
# MATLAB: a = 0.42748*R^2*Tc^2.5/Pc; b = 0.08664*R*Tc/Pc;
a = 0.42748 * (R**2) * (Tc**2.5) / Pc
b = 0.08664 * R * Tc / Pc

# 3. SRK 상태 방정식 정의
# f(V) = R*T/(V-b) - a/(V*(V+b)*sqrt(T)) - P = 0
def srk_eqn(V):
    term1 = (R * T) / (V - b)
    term2 = a / (V * (V + b) * np.sqrt(T))
    return term1 - term2 - P

# 4. 해 구하기
V0 = 0.1  # 초기 추측값 (Initial guess)
V_solution = fsolve(srk_eqn, V0)

print(f"SRK 비체적 (Specific Volume, V): {V_solution[0]:.6f} L/mol")