import numpy as np

# 데이터 설정 (MATLAB 원본 상수)
d = 0.614      # 직경 (inch)
Lst = 20       # 길이 (ft)
f = 0.026      # 마찰 계수
K = 2.026      # 부차적 손실 계수
Z = 0.9        # 압축 인자
T_c = 100      # 온도 (F)
Mw = 19.5      # 분자량
P1 = 1124.7    # 입구 압력 (psia)
P2 = 414.7     # 출구 압력 (psia)
k = 1.27       # 비열비
mu = 0.012     # 점도 (lb/ft/hr 단위로 추정됨)

# 단위 변환 및 보조 계산
T = T_c + 460  # 화씨를 랭킨(Rankine)으로 변환
D = d / 12     # 직경을 ft로 변환
Area = np.pi * D**2 / 4

# 밀도 계산 (lb/ft^3)
rho = (P1 * Mw) / (10.72 * Z * T)

# 압력 강하 및 손실 계수
delP = P1 - P2
Kp = f * Lst / D
Kt = K + Kp    # 총 손실 계수

# 질량 유량 (G) 계산 (lb/hr 단위)
# 공식: G = 1335.6 * d^2 * sqrt(rho * (P1^2 - P2^2) / (Kt + 2 * ln(P1/P2)) / P1)
G = 1335.6 * d**2 * np.sqrt(rho * (P1**2 - P2**2) / (Kt + 2 * np.log(P1 / P2)) / P1)

# 유속, 마하수 및 임계 압력 계산
Nre = 6.31 * G / (D * mu)            # 레이놀즈 수
Vg = 0.0509 * G / (rho * d**2)       # 유속 (ft/s)
Vs = 223 * np.sqrt(k * T / Mw)       # 음속 (ft/s)
Mach = Vg / Vs                       # 마하수
Sg = Mw / 29                         # 비중
R = 1544 / (29 * Sg)                 # 기체 상수
Pc = (G / (11400 * d**2)) * np.sqrt(R * T / (k * (k + 1))) # 임계 압력

# 결과 출력
print(f"G (Mass Flow): {G:.4f} lb/hr")
print(f"Vg (Velocity): {Vg:.4f} ft/s")
print(f"Mach Number: {Mach:.4f}")
print(f"Pc (Critical Pressure): {Pc:.4f}")