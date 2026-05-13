import numpy as np

# 1. 변수 정의 (영국 단위계 기반)
rho = 62.4      # 밀도 (lb/ft^3)
mu = 6.905e-4   # 점도 (lb/ft·s)
D = 1/3         # 직경 (ft)
Q = 0.2         # 유량 (ft^3/s)
K = 0.45 + 3*0.5 + 1  # 부차적 손실 계수
gc = 32.174     # 중력 변환 상수
L = 5 + 300 + 100 + 120 + 20 # 총 배관 길이 (ft)
eD = 5e-5       # 상대 조도 (epsilon/D)

# 2. 기초 물성 계산
v = 4 * Q / (np.pi * D**2)      # 유속 (ft/s)
m = Q * rho                     # 질량 유량 (lb/s)
Nre = D * v * rho / mu          # 레이놀즈 수

# 3. 마찰 계수(f) 계산 (조건문 사용)
if Nre < 2100:
    f = 16.0 / Nre
else:
    # Shacham 방정식 적용
    term1 = eD / 3.7
    term2 = -5.02 * np.log10(eD / 3.7 + 14.5 / Nre) / Nre
    den = 16 * (np.log10(term1 + term2))**2
    f = 1.0 / den

# 4. 에너지 손실 및 일률(Work) 계산
# Darcy-Weisbach 식 및 에너지 평형식
dH = (4 * f * L / D + K) * (v**2) / (2 * gc)
Ws = m * (dH + (105 - 20))

# 5. 결과 출력
print(f"Velocity (v): {v:.4f} ft/s")
print(f"Friction factor (f): {f:.6f}")
print(f"Work/Power (Ws): {Ws:.4f} ft·lb/s")