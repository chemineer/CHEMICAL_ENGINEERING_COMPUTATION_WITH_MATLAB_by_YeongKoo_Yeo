import numpy as np
from scipy.optimize import fsolve

# 1. 입력 데이터 설정
e_D = 1.3e-4    # 상대 거칠기 (Relative roughness)
Nre = 6.5e4     # 레이놀즈 수 (Reynolds number)
f0 = 0.1        # 초기 추측값 (Initial guess)

# 2. Colebrook 방정식 정의
# 1/sqrt(f) + 0.86 * ln(eD/3.7 + 2.51/(Nre * sqrt(f))) = 0 형태
def colebrook(f):
    # f가 0 이하가 되지 않도록 안전장치를 두거나 sqrt 내부를 보호합니다.
    term1 = 1 / np.sqrt(f)
    term2 = 0.86 * np.log(e_D / 3.7 + 2.51 / (Nre * np.sqrt(f)))
    return term1 + term2

# 3. fsolve를 사용하여 해 구하기
f_solution = fsolve(colebrook, f0)

print(f"마찰 계수 (Friction Factor, f): {f_solution[0]:.6f}")