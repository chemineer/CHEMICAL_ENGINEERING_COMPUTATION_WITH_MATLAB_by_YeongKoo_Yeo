import numpy as np

# 데이터 설정
G = 582 / 3600    # 질량 유량 (kg/s)
d1 = 0.10226      # 상류 배관 직경 (m)
d2 = 0.01         # 오리피스 직경 (m)
p1 = 1.2e6        # 상류 압력 (Pa)
Fa = 1            # 열팽창 계수
rho = 10.25       # 밀도 (kg/m^3)
mu = 1.3e-5       # 점도 (Pa*s)
gam = 1.3         # 비열비
K = 1             # 손실 계수

# 기본 계산
beta = d2 / d1
A1 = np.pi * d1**2 / 4
q1 = G / rho
v1 = q1 / A1      # 상류 속도 (m/s)
Ao = np.pi * d2**2 / 4  # 오리피스 면적
Nre = d1 * v1 * rho / mu  # 레이놀즈 수
Cd = 0.73         # 임계 흐름 유량 계수

# 팽창 계수 (Y) 계산
# 공식: G / (Cd * Fa * Ao) / sqrt(2 * rho * p1 * gam * (2/(1+gam))^((gam+1)/(gam-1)))
denom = np.sqrt(2 * rho * p1 * gam * (2 / (1 + gam))**((gam + 1) / (gam - 1)))
Y = G / (Cd * Fa * Ao) / denom

# 압력 강하 (dp) 계산
# 공식: p1 / (0.9953 + 0.9054/sqrt(K) + 0.1173/K - 0.0195/K^1.5)
dp = p1 / (0.9953 + 0.9054 / np.sqrt(K) + 0.1173 / K - 0.0195 / K**1.5)

# 결과 출력
print(f"Expansion factor = {Y:.6g}")
print(f"Pressure drop = {dp / 1000:.6g} kPa")