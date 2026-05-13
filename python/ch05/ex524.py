import numpy as np

# 데이터 설정
T = 60          # 온도 (C)
sg = 0.42       # 비중
d = 0.0779      # 직경 (m)
K = 2.848       # 총 손실 계수
p0 = 101.3      # 기준 압력 (kPa)
p1 = p0 + np.array([800, 150])  # 실린더 압력 (kPa)
n = len(p1)

# 한계 압력 강하 비 (limiting dp/p1) 계산
dp1c = 1 / (0.9953 + 0.9054 / np.sqrt(K) + 0.1173 / K - 0.0195 / K**1.5)

# 추정 압력 강하 비
dp1e = (p1 - p0) / p1

# 한계 팽창 계수 (limiting expansion factor)
Y_limit = 0.0415 * np.log(K) + 0.6097

# 루프를 통한 유량(Q) 계산
for k in range(n):
    # 흐름이 음속(sonic)인 경우
    if dp1c < dp1e[k]:
        Y = Y_limit
        dp = dp1c * p1[k]
        # 유량 계산 (m^3/hr)
        Q = 3600 * 53.64 * Y * d**2 * np.sqrt(dp * p1[k] / (K * (T + 273.15) * sg))
        print(f"Flow rate (cylinder pressure: {p1[k]-p0:.1f} kPaG) = {Q:.4f} m^3/hr")
        
    # 흐름이 아음속(subsonic)인 경우
    else:
        m = (1 - Y_limit) / dp1c  # 기울기
        Y = 1 - m * dp1e[k]       # 비질식 유동(unchoked flow)
        dp = (p1[k] - p0)         # 실제 압력 강하
        Q = 3600 * 53.64 * Y * d**2 * np.sqrt(dp * p1[k] / (K * (T + 273.15) * sg))
        print(f"Flow rate (cylinder pressure: {p1[k]-p0:.1f} kPaG) = {Q:.4f} m^3/hr")