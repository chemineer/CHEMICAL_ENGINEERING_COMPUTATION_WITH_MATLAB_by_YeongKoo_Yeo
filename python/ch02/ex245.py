import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import trapezoid

# 1. 데이터 정의
t = np.array([20.1, 49.8, 81.5, 109.8, 140.2, 171.1, 201.4, 229.5])  # 온도 (deg.C)
Cp = np.array([28.98, 29.11, 29.96, 29.47, 29.67, 29.91, 30.02, 30.14])  # Cp (J/mol/deg.C)

n = 6.5       # 몰 수 (mol)
t1 = 55       # 시작 온도
t2 = 185      # 종료 온도

# --- Step 1: 3차 스플라인 보간을 사용하여 Cp1, Cp2 찾기 ---
# MATLAB의 interp1(..., 'spline')과 동일한 역할을 수행합니다.
cs = CubicSpline(t, Cp)
Cp1 = cs(t1)
Cp2 = cs(t2)

# --- Step 2: 수치 적분을 위한 하위 구간 데이터 구성 ---
# t1과 t2 사이의 원본 데이터 포인트들을 추출하고 보간된 값을 양 끝에 추가합니다.
ts = np.array([t1, 81.5, 109.8, 140.2, 171.1, t2])
Cps = np.array([Cp1, 29.96, 29.47, 29.67, 29.91, Cp2])

# --- Step 3: 사다리꼴 공식을 이용한 수치 적분 수행 ---
# MATLAB의 trapz(ts, Cps)를 수행합니다. (y, x 순서 주의)
integral_value = trapezoid(Cps, ts)
delH = n * integral_value

# 결과 출력
print(f"Cp at t1 ({t1}): {Cp1:.4f}")
print(f"Cp at t2 ({t2}): {Cp2:.4f}")
print(f"Delta H (Enthalpy Change): {delH:.4f} J")