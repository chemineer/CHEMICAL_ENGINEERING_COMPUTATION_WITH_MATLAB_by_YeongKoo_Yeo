import numpy as np
from scipy.optimize import root_scalar

# 1. 상수 정의
rho = 997.92
Q = 0.0085
mu = 0.000982
L = 250
gc = 1
Ff = 27.5

# 2. 방정식 정의 (g(D) = 0)
def g(D):
    # D가 0 이하가 되는 경우를 방지 (물리적 범위 제한)
    if D <= 1e-6:
        return 1e6 
    
    # 레이놀즈 수 및 속도 수두 계산
    reynolds = (4 * rho * Q) / (np.pi * mu * D)
    velocity_head = (4 * Q / (np.pi * D**2))**2
    
    # 방정식 값 반환
    return 0.0936 * (reynolds**(-0.2)) * (L / D / gc) * velocity_head - Ff

# 3. 개선된 수치 해석 (브래킷 방식 사용)
# fsolve doesn't work properly
# 물리적으로 가능한 직경 범위를 1mm(0.001m)에서 5m 사이로 설정
bracket = [0.001, 5.0]

try:
    sol = root_scalar(g, bracket=bracket, method='brentq')
    
    if sol.converged:
        print(f"해를 찾았습니다: D = {sol.root:.6f} m")
    else:
        print("수렴에 실패했습니다. 구간 범위를 조정해 보세요.")
        
except ValueError as e:
    print(f"오류 발생: 구간 내에 해가 존재하지 않을 수 있습니다. ({e})")