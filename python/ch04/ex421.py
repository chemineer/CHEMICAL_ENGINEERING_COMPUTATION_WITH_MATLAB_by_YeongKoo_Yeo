import numpy as np
from scipy.optimize import fsolve

# --- 데이터 설정 ---
zf = np.array([0.1, 0.25, 0.5, 0.15])
A = np.array([6.64380, 6.82915, 6.80338, 6.80776])
B = np.array([395.74, 663.72, 804.00, 935.77])
C = np.array([266.681, 256.681, 247.040, 238.789])
P_mmHg = np.array([16, 18, 20, 24]) * 760

# --- 함수 정의 (스칼라 값 처리 추가) ---

def sumf(alp, zf, k):
    # alp[0]으로 명시하여 첫 번째 값을 추출
    a = alp[0] if hasattr(alp, "__len__") else alp
    res = np.sum(zf * (1 - k) / (1 + a * (k - 1)))
    return res

def zdk(T, zf, P, A, B, C):
    T_val = T[0] if hasattr(T, "__len__") else T
    Pv = 10**(A - B / (C + T_val))
    k = Pv / P
    return np.sum(zf / k) - 1

def zmk(T, zf, P, A, B, C):
    T_val = T[0] if hasattr(T, "__len__") else T
    Pv = 10**(A - B / (C + T_val))
    k = Pv / P
    return np.sum(zf * k) - 1

# --- 메인 계산 ---
T_fixed = 50.0
Pv = 10**(A - B / (C + T_fixed))

print(f"{'P (atm)':<10} | {'Alpha':<10} | {'x_1':<8} | {'y_1':<8}")
print("-" * 50)

for p in P_mmHg:
    k = Pv / p
    # fsolve 결과의 첫 번째 요소[0]를 사용하여 스칼라 추출
    ax = fsolve(sumf, 0.5, args=(zf, k))[0]
    xv = zf / (1 + ax * (k - 1))
    yv = k * xv
    print(f"{p/760:<10.1f} | {ax:<10.4f} | {xv[0]:<8.4f} | {yv[0]:.4f}")

# Tdew, Tbub 계산
Tdew = [fsolve(zdk, 50.0, args=(zf, p, A, B, C))[0] for p in P_mmHg]
Tbub = [fsolve(zmk, 50.0, args=(zf, p, A, B, C))[0] for p in P_mmHg]

print("\n--- 결과 ---")
print(f"Dew Point Temperatures: {np.round(Tdew, 2)}")
print(f"Bubble Point Temperatures: {np.round(Tbub, 2)}")