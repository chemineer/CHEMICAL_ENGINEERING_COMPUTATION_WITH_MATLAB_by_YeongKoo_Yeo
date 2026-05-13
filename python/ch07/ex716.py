import math

# 1. 데이터 입력 (Data)
Ht = 0.6      # 전달 단위 높이 (Height of a transfer unit)
H = 26        # 헨리 상수 (Henry's law constant)
Gm = 206      # 가스 몰 유량 (Gas molar flow rate)
Lm = 12240    # 액체 몰 유량 (Liquid molar flow rate)
x2 = 0.0      # 입구 액체 내 SO2 농도
y1 = 0.03     # 입구 가스 내 SO2 농도
y2 = 0.003    # 출구 가스 내 SO2 농도

# 2. 계산 (Calculations)
# 흡수 인자 (Absorption factor)
Sf = Lm / (H * Gm)

# 이론적 전달 단위 수 (Number of transfer units, Nt)
# MATLAB의 log()는 자연로그이므로 math.log()를 사용합니다.
term1 = Sf / (Sf - 1)
term2 = (1 - 1/Sf) * ((y1 - H*x2) / (y2 - H*x2)) + 1/Sf
Nt = term1 * math.log(term2)

# 전체 충전 높이 (Total packing height, Hpd)
Hpd = Nt * Ht

# 3. 결과 출력
print(f"Number of theoretical transfer units = {Nt:g}")
print(f"Total packing height = {Hpd:g} m")