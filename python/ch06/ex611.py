import numpy as np
from scipy.optimize import fsolve

# 1. 데이터 및 파라미터 설정
F0 = 40        # 유입 유량
F = 40         # 유출 유량
Fj = 49.9      # 자켓 냉각수 유량
Ca0 = 0.55     # 유입 농도
V = 48         # 반응기 부피
rho = 50       # 반응액 밀도
rhoj = 62.3    # 냉각수 밀도
Cp = 0.75      # 반응액 비열
Cj = 1         # 냉각수 비열
A = 250        # 열전달 면적
U = 150        # 총괄 열전달 계수
T0 = 530       # 유입 온도 (Rankine 또는 Kelvin)
Tj0 = 530      # 자켓 유입 온도
alp = 7.08e10  # 빈도 계수 (Pre-exponential factor)
lam = -3e4     # 반응열 (Heat of reaction)
E = 3e4        # 활성화 에너지
R = 1.9872     # 기체 상수

# 2. 비선형 방정식 시스템 정의
def exocstr(x):
    # x[0] = Ca, x[1] = T, x[2] = Tj
    Ca, T, Tj = x
    
    # 아레니우스 속도 항 계산
    rate_term = alp * V * Ca * np.exp(-E / (R * T))
    
    # f1: 물질 수지 (Mass Balance)
    f1 = F0 * Ca0 - F * Ca - rate_term
    
    # f2: 에너지 수지 (Energy Balance - Reactor)
    f2 = rho * Cp * (F0 * T0 - F * T) - lam * rate_term - U * A * (T - Tj)
    
    # f3: 에너지 수지 (Energy Balance - Jacket)
    f3 = rhoj * Cj * Fj * (Tj0 - Tj) + U * A * (T - Tj)
    
    return [f1, f2, f3]

# 3. 초기 추정값 설정 및 해 계산
x0 = [Ca0, T0, Tj0]  # 초기 추정값: 유입 조건과 동일하게 설정
solution = fsolve(exocstr, x0)

# 4. 결과 출력
Ca_ss, T_ss, Tj_ss = solution
print(f"--- 정상 상태 결과 ---")
print(f"농도 (Ca): {Ca_ss:.4f}")
print(f"반응기 온도 (T): {T_ss:.2f}")
print(f"자켓 온도 (Tj): {Tj_ss:.2f}")