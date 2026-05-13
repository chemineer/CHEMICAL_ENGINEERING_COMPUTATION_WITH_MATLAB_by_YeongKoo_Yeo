import numpy as np
from scipy.optimize import fsolve

# 1. 데이터 설정
dP = -18 * 1.01325e5
dz = 100
rho = 1e3
mu = 1e-3
L = 1500
D = 0.154
rh = 4.57e-5
g = 9.8

# 2. 연립 방정식 정의
def qdpfun(x, dP, dz, rho, mu, L, D, rh, g):
    v = x[0]
    f = x[1]
    
    # 0으로 나누기 방지
    if v <= 0 or f <= 0:
        return [1e6, 1e6]
    
    Nre = D * v * rho / mu  # 레이놀즈 수
    
    # 방정식 1: 에너지 평형식
    # -v^2/2 + g*dz + dP/rho + 2*f*L*v^2/D = 0
    eq1 = -v**2 / 2 + g * dz + dP / rho + 2 * f * L * v**2 / D
    
    # 방정식 2: 콜브룩-화이트 식 (1/sqrt(f) = -2*log10(rh/3.7/D + 2.51/Nre/sqrt(f)))
    # 제공된 코드 식: 1/sqrt(f) + 1.7372*log(rh/3.7/D + 1.255/Nre/sqrt(f)) = 0
    # 참고: 파이썬의 np.log는 자연로그(ln)입니다.
    eq2 = 1 / np.sqrt(f) + 1.7372 * np.log(rh / (3.7 * D) + 1.255 / (Nre * np.sqrt(f)))
    
    return [eq1, eq2]

# 3. 비선형 방정식 풀이
x0 = [5, 0.001]  # 초기값
sol = fsolve(qdpfun, x0, args=(dP, dz, rho, mu, L, D, rh, g))

v, f = sol[0], sol[1]
Q = v * (np.pi * D**2) / 4  # 유량
dPf = 2 * f * rho * L * v**2 / D  # 마찰 손실 압력 강하

# 4. 결과 출력
print(f"Volumetric flow rate of water = {Q:.6g} m^3/sec")
print(f"Pressure drop due to friction loss = {dPf/1000:.6g} kPa")