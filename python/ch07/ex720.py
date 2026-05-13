import numpy as np
from scipy.optimize import fsolve

# 1. 데이터 입력 (Feed composition & Antoine constants)
# 1: ethane, 2: propane(LK), 3: i-butane(HK), 4: n-butane
xF = np.array([0.15, 0.18, 0.18, 0.49])
A = np.array([3.93835, 4.53678, 4.3281, 4.35576])
B = np.array([659.739, 1149.36, 1132.108, 1175.581])
C = np.array([-16.719, 24.906, 0.918, -2.071])
P = 7

# 2. Distillate(D) 및 Bottom(B) 조성 계산
Br = np.zeros(4)
Br[0] = 0
Br[2] = 0.95 * xF[2]
Br[3] = xF[3]
Br[1] = 0.001 * (Br[2] + Br[3]) / 0.999

Dr = np.zeros(4)
Dr[0] = xF[0]
Dr[1] = xF[1] - Br[1]
Dr[2] = 0.05 * xF[2]
Dr[3] = 0

xD = Dr / np.sum(Dr)
xB = Br / np.sum(Br)

# 3. 이슬점(D) 및 버블점(B) 계산
Td0 = 270
Tb0 = 300

# fD: Td에 대한 비선형 방정식 (이슬점)
def fD_func(Td):
    Psat = 10**(A - B / (Td + C))
    return np.sum(P * xD / Psat) - 1

# fB: Tb에 대한 비선형 방정식 (버블점)
def fB_func(Tb):
    Psat = 10**(A - B / (Tb + C))
    return np.sum(xB * Psat / P) - 1

Td = fsolve(fD_func, Td0)[0]
Tb = fsolve(fB_func, Tb0)[0]

# 4. Fenske 식 (최소 단수 Nm 계산)
kD = (10**(A - B / (Td + C))) / P
kB = (10**(A - B / (Tb + C))) / P

num = np.log((xD[1] / xD[2]) * (xB[2] / xB[1]))
den = np.log(np.sqrt((kD[1] / kD[2]) * (kB[1] / kB[2])))
Nm = num / den

# 5. 원료 공급단 온도(Tf) 계산
Tf0 = 300
fF_func = lambda Tf: np.sum(xF * 10**(A - B / (Tf + C)) / P) - 1
Tf = fsolve(fF_func, Tf0)[0]
kF = (10**(A - B / (Tf + C))) / P

# 6. Underwood 식 (최소 환류비 Rm 계산)
th0 = 1.5
def fTh_func(theta):
    # kF[2]는 i-butane(HK)의 k값
    return np.sum(xF * kF / (kF - theta * kF[2]))

theta = fsolve(fTh_func, th0)[0]
alpD = kD / kD[2]
Rm = np.sum(alpD * xD / (alpD - theta)) - 1

# 7. Gilliland 식 (실제 단수 N 계산)
R = 1.5 * Rm
X_gill = (R - Rm) / (R + 1)  # MATLAB 코드의 (1.5-1)*Rm/(R+1)와 동일
Y_gill = 0.75 * (1 - X_gill**0.5658)
N = (Y_gill + Nm) / (1 - Y_gill)

# 8. Kirkbride 식 (원료 공급단 위치 계산)
c1 = np.sum(Br) / np.sum(Dr)
c2 = xF[2] / xF[1]
c3 = (xB[1] / xD[2])**2
h = (c1 * c2 * c3)**0.206
p = N / (1 + h)  # 공급단 아래 단수
m = p * h       # 공급단 위 단수

# 9. 결과 출력
print(f"Bubble point(B) = {Tb:8.4f}, Dew point(D) = {Td:8.4f}")
print(f"Feed temp.(F) = {Tf:8.4f}")
print(f"theta = {theta:8.4f}, min. number of stages = {Nm:7.4f}")
print(f"Min. reflux ratio = {Rm:7.4f}, actual number of stages = {N:7.4f}") # MATLAB 출력 순서 조정
print(f"Stages above the feed stage = {m:7.4f}")
print(f"Stages below the feed stage = {p:7.4f}")