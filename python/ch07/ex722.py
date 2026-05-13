import numpy as np
from scipy.optimize import fsolve

# 1. 데이터 설정
# 성분: 1: n-hexane, 2: n-heptane, 3: n-octane
xf = np.array([0.33, 0.37, 0.3])    # 피드 몰 분율
xd = np.array([0.99, 0.01, 0.0])   # 유출액 몰 분율
xb_1 = 0.01                        # 하부 제품 내 n-hexane 분율
F = 100                            # 피드 유량 (mol/h)
Tf = 105                           # 피드 온도 (deg.C)
P = 1.2                            # 압력 (atm)
q = 0.4                            # 피드 상태 (60% 증기)

# Antoine 상수 (A, B, C)
A = np.array([6.87024, 6.89385, 6.90940])
B = np.array([1168.720, 1264.370, 1349.820])
C = np.array([224.210, 216.636, 209.385])

# 2. 물질 수지 계산 (Product flow rates and composition)
# z[0]=Distillate 유량, z[1]=xb(2), z[2]=xb(3)
def material_balance(z):
    Dist = z[0]
    xb = np.array([xb_1, z[1], z[2]])
    Bott = F - Dist
    # 성분별 수지: F*xf_i = D*xd_i + B*xb_i
    f = F * xf - Dist * xd - Bott * xb
    return f

z0 = [50, 0.5, 0.5]
z_sol = fsolve(material_balance, z0)
Dist = z_sol[0]
xb = np.array([xb_1, z_sol[1], z_sol[2]])
Bott = F - Dist

# 3. 끓는점 및 이슬점 계산
# 유출액(D)의 이슬점(Dew point) Td
def dew_point_eq(Td):
    psat = 10**(A - B / (Td + C))
    return (760 * P) * np.sum(xd / psat) - 1

# 하부액(B)의 버블점(Bubble point) Tb
def bubble_point_eq(Tb):
    psat = 10**(A - B / (Tb + C))
    return np.sum(xb * psat) / (760 * P) - 1

Td = fsolve(dew_point_eq, 100)[0]
Tb = fsolve(bubble_point_eq, 100)[0]

# 4. 평형 계수(K) 및 상대 휘발도(alpha) 계산
def get_K(T):
    return 10**(A - B / (T + C)) / (760 * P)

Kf = get_K(Tf)
Kd = get_K(Td)
Kb = get_K(Tb)

# 상대 휘발도 (LK/HK 기준, 여기서는 1번과 2번 성분)
alpaf12 = Kf[0] / Kf[1]
alpad12 = Kd[0] / Kd[1]
alpab12 = Kb[0] / Kb[1]
avgalpa = (alpaf12 * alpad12 * alpab12)**(1/3)

# 5. 결과 계산
# (1) 최소 단수 (Fenske eqn.)
xlkd, xhkd = xd[0], xd[1]
xlkb, xhkb = xb[0], xb[1]
Nmin = np.log((xlkd / xhkd) * (xhkb / xlkb)) / np.log(avgalpa)

# 최소 환류비 (Underwood eqn.)
alpaf = Kf / Kf[1]
alpad = Kd / Kd[1]

def underwood_theta(theta):
    return 1 - q - np.sum(alpaf * xf / (alpaf - theta))

theta_sol = fsolve(underwood_theta, 1.5)[0]
Rmin = np.sum(alpad * xd / (alpad - theta_sol)) - 1
Ract = 1.5 * Rmin

# (2) 이론 단수 (Gilliland/Eduljee correlation)
X = (Ract - Rmin) / (Ract + 1)
Y = 0.75 * (1 - X**0.5668)
N = np.ceil((Nmin + Y) / (1 - Y))

# (3) 피드 단 위치 (Kirkbride eqn.)
xlkf, xhkf = xf[0], xf[1]
Nrs = (xhkf * xlkb**2 * Bott / (xlkf * xhkd**2 * Dist))**0.206
Ns = np.ceil(N / (Nrs + 1))
Nr = N - Ns
Ftray = Nr

# 결과 출력
print(f"Dew point(D) = {Td:.4f} (deg.C), bubble point(B) = {Tb:.4f} (deg.C)")
print(f"(1) the minimum number of trays = {Nmin:.4f}")
print(f"(2) the number of ideal trays (R=1.5*Rmin) = {N}")
print(f"(3) feed tray location = {Ftray}")