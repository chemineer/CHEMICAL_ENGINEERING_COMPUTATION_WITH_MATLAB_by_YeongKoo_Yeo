import numpy as np
from scipy.optimize import fsolve

# 1. 시스템 방정식 정의 (변화율이 0이 되는 지점을 찾음)
def cstrst(x, Ca0, k0, E, tau, Tc, kappa, Cp):
    # x[0] = Ca (농도), x[1] = T (온도)
    Ca = x[0]
    T = x[1]
    
    # 반응 속도 ra 및 반응열 dHr 계산
    ra = -k0 * np.exp(-E / T) * Ca
    dHr = -151080 + 2 * (T - 298.15)
    
    # 두 방정식의 값이 0이 되어야 함
    # f1: 농도 수지, f2: 에너지 수지
    f1 = (Ca0 - Ca) / tau - (-ra)
    f2 = (-dHr / Cp) * (-ra / Ca0) - (1 + kappa) * (T - Tc) / tau
    
    return [f1, f2]

# 2. 데이터 및 파라미터 설정
Ca0 = 0.4
k0 = 460
E = 1380
tau = 0.18
Tc = 298.15
kappa = 78
Cp = 32

# 3. 초기 추정값 및 해 구하기
x0 = [0.1, 300]  # 초기 추정값 [Ca, T]
solution = fsolve(cstrst, x0, args=(Ca0, k0, E, tau, Tc, kappa, Cp))

# 4. 결과 출력 (K -> deg.C 변환 포함)
Ca_ss = solution[0]
T_ss_C = solution[1] - 273.15

print(f"정상 상태 농도 (Ca): {Ca_ss:.4f} mol/cm^3")
print(f"정상 상태 온도 (T): {T_ss_C:.2f} °C")