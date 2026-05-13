import numpy as np

# 1. 데이터 입력 (Data)
rhog = 1.17      # 가스 밀도
rhol = 1000      # 액체 밀도
Fp = 130         # 충전물 계수
phi = 1          # 액체 밀도 보정 계수
mul = 0.8        # 액체 점도
Lm = 2450        # 최소 액체 유량
G = 103          # 가스 유량
g = 9.82         # 중력 가속도
f = 0.7          # 플러딩(Flooding) 대비 실제 속도 비율

# 2. 직경 계산 (Find diameter)
L = 1.5 * Lm
X = (L / G) * np.sqrt(rhog / rhol)

# 로그 계산 (MATLAB의 log는 자연로그인 np.log와 동일)
log_X = np.log(X)
z = -1.668 - 1.085 * log_X - 0.297 * (log_X**2)
Y = 10**z

# 범람 속도(Superficial velocity at flooding) 및 실제 속도 계산
Gf = np.sqrt(rhog * rhol * g * Y / (phi * Fp * mul**0.2))
Gr = 0.7 * Gf

# 단면적 및 직경 계산
A = (G / 60) / Gr
Dt = np.sqrt(4 * A / np.pi)

# 3. 결과 출력
print(f"Column diameter = {Dt:g} m")