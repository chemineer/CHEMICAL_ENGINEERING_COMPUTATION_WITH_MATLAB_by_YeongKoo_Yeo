import math

# 데이터 설정 (MATLAB의 cfdat 배열)
# cfdat는 14행 3열의 0으로 채워진 배열로 시작함
cfdat = [[0.0 for _ in range(3)] for _ in range(14)]
cfdat[1][0] = 6.0
cfdat[10][0] = 2.0

col2_data = [10, 8, 5, 3, 5, 10, 15, 10, 8, 15, 1.5, 1, 1.6, 0]
col3_data = [3, 2, 1.5, 1, 1.5, 2.5, 40, 20, 2.5, 15, 0.5, 0, 5, 10]

for i in range(14):
    cfdat[i][1] = 100 * col2_data[i]
    cfdat[i][2] = 0.1 * col3_data[i]

# 상수 및 기본 변수 설정
g = 32.2
d = 6.065
D = d / 12
pr = 1.5e-4
w = 75000
mu = 1.25
rho = 64.3

Lst = 78
z = 8
Area = math.pi * D**2 / 4
Nre = 6.31 * w / (d * mu)
v = 0.0509 * w / (rho * d**2)
Q = w / (8.02 * rho)

# 마찰 계수 (f) 계산
if Nre <= 2100:
    f = 64 / Nre
else:
    # Chen 방정식을 이용한 난류 마찰 계수
    Av = pr / (3.7 * D) + (6.7 / Nre)**0.9
    f = 4.0 / (-4 * math.log10(pr / D / 3.7 - 5.02 * math.log10(Av) / Nre))**2

# Ksum, Klsum, Kt 계산
Ksum = sum(cfdat[i][0] * cfdat[i][1] for i in range(12))
Klsum = sum(cfdat[i][0] * cfdat[i][2] for i in range(12))
Kt = Ksum / Nre + Klsum * (1 + 1 / d)

# Kee1, Kee2, K, Leq, L 계산
Kee1 = sum(cfdat[i][0] * cfdat[i][1] for i in range(12, 14))
Kee2 = sum(cfdat[i][0] * cfdat[i][2] for i in range(12, 14))
K = Kee1 / Nre + Kee2 + Kt
Kt = Kt + K + f * Lst / D
Leq = K * D / f
L = Lst + Leq

# 압력 강하(delP) 및 손실(delPsi, delH) 계산
if Nre <= 2100:
    delP = 0.0034 * mu * w / (d**4 * rho)
else:
    delP = 0.000336 * f * w**2 / (d**5 * rho)

delPsi = delP * L / 100 + z * rho / 144
delH = 0.000483 * f * L * w**2 / (d**5 * rho**2) + z

# 결과 출력
print(f"Leq: {Leq}")
print(f"L: {L}")
print(f"delPsi: {delPsi}")
print(f"delH: {delH}")