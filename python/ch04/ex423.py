import numpy as np
from unifgam import unifgam

# 1. 데이터 설정
n = 4
P = 1.0  # atm
Rg = 82.06  # cm^3 atm/(mol K)
x = np.array([0.162, 0.068, 0.656, 0.114])

# Antoine 상수
A = np.array([9.2033, 12.2786, 9.1690, 9.2675])
B = np.array([2697.55, 3803.98, 2731.00, 2788.51])
C = np.array([-48.78, -41.68, -47.11, -52.36])

# 비리알 계수 행렬 (대칭 행렬 구성)
Bij = np.zeros((n, n))
# MATLAB 코드의 typo (Bij(1,2) -657.0)를 고려하여 할당
Bij[0,0] = -1360.1; Bij[0,1] = -657.0; Bij[0,2] = -1274.2; Bij[0,3] = -1218.8
Bij[1,1] = -1174.7; Bij[1,2] = -621.8; Bij[1,3] = -589.7
Bij[2,2] = -1191.9; Bij[2,3] = -1137.9
Bij[3,3] = -1086.9

# 대칭화 (MATLAB: for i = 1:n-1, Bij(:,i) = Bij(i,:); end)
for i in range(n):
    for j in range(i + 1, n):
        Bij[j, i] = Bij[i, j]

# UNIFAC용 상수
k, R, Q = 6, np.array([0.9011, 0.6744, 0.4469, 0.2195, 0.5313, 1.0000]), np.array([0.848, 0.540, 0.228, 0.000, 0.400, 1.200])
nu = np.array([[2, 1, 3, 0], [4, 1, 1, 0], [0, 0, 1, 0], [0, 0, 1, 0], [0, 0, 0, 6], [0, 1, 0, 0]])
amn = np.array([[0, 0, 0, 0, 61.13, 986.5], [0, 0, 0, 0, 61.13, 986.5], [0, 0, 0, 0, 61.13, 986.5],
                [0, 0, 0, 0, 61.13, 986.5], [-11.12, -11.12, -11.12, -11.12, 0, 636.1],
                [156.4, 156.4, 156.4, 156.4, 89.60, 0]])

# Delta 행렬 계산
delta = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        delta[j, i] = 2 * Bij[j, i] - Bij[i, i] - Bij[j, j]
        if j == i: delta[i, j] = 0

# Step 1) 초기 온도 추정
PHIi = np.ones(n)
Tsat = B / (A - np.log(P)) - C
T = np.sum(x * Tsat)

# Step 2) 초기 루프 전 계산 (MATLAB 로직과 동기화)
Psat = np.exp(A - B / (T + C))
gamma = unifgam(k, R, Q, nu, amn, n, x, T)
P1sat = P / np.sum(x * gamma * Psat / (PHIi * Psat[0]))
T = B[0] / (A[0] - np.log(P1sat)) - C[0]

# Step 3) 반복 계산
Told = T
Terr, Tcrit = 1.0, 1e-6

while Terr > Tcrit:
    Psat = np.exp(A - B / (Told + C))
    # 기상 조성 y 계산
    y = (x * gamma * Psat) / (PHIi * P)
    
    # PHI 계산 (Triple loop)
    for i in range(n):
        sumy = 0.0
        for j in range(n):
            for ki in range(n):
                sumy += y[j] * y[ki] * (2 * delta[j, i] - delta[j, ki])
        PHIi[i] = np.exp((Bij[i, i] * (P - Psat[i]) + sumy * P / 2) / (Rg * Told))
    
    # gamma 갱신 및 새로운 T 계산
    gamma = unifgam(k, R, Q, nu, amn, n, x, Told)
    P1sat = P / np.sum(x * gamma * Psat / (PHIi * Psat[0]))
    T = B[0] / (A[0] - np.log(P1sat)) - C[0]
    
    Terr = abs(Told - T)
    Told = T

# 결과 출력
print(f"최종 온도 T: {T:.4f}")
print(f"기상 조성 y: {y}")
print(f"포화 증기압 Psat: {Psat}")
print(f"PHIi: {PHIi}")
print(f"gamma: {gamma}")