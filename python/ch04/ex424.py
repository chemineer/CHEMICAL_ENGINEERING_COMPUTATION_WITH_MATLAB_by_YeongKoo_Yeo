import numpy as np
from scipy.optimize import fsolve
from unifgam import unifgam

# 1. 데이터 설정
n = 4
P = 1.0  # atm
T = 334.15  # K
Rg = 82.06
z = np.array([0.25, 0.40, 0.20, 0.15])

A = np.array([9.2033, 12.2786, 9.1690, 9.2675])
B = np.array([2697.55, 3803.98, 2731.00, 2788.51])
C = np.array([-48.78, -41.68, -47.11, -52.36])

# 비리알 계수 행렬 (0-indexing, 대칭화 처리 필요)
Bij = np.zeros((n, n))
# 데이터 할당 (MATLAB 오타 수정)
Bij[0,0] = -1360.1; Bij[0,1] = -657.0; Bij[0,2] = -1274.2; Bij[0,3] = -1218.8
Bij[1,1] = -1174.7; Bij[1,2] = -621.8; Bij[1,3] = -589.7
Bij[2,2] = -1191.9; Bij[2,3] = -1137.9; Bij[3,3] = -1086.9
for i in range(n):
    for j in range(i+1, n):
        Bij[j, i] = Bij[i, j]

# UNIFAC용 파라미터 (이전과 동일)
k, R, Q = 6, np.array([0.9011, 0.6744, 0.4469, 0.2195, 0.5313, 1.0000]), np.array([0.848, 0.540, 0.228, 0.000, 0.400, 1.200])
nu = np.array([[2,1,3,0], [4,1,1,0], [0,0,1,0], [0,0,1,0], [0,0,0,6], [0,1,0,0]])
amn = np.array([[0,0,0,0,61.13,986.5], [0,0,0,0,61.13,986.5], [0,0,0,0,61.13,986.5],
                [0,0,0,0,61.13,986.5], [-11.12,-11.12,-11.12,-11.12,0,636.1], [156.4,156.4,156.4,156.4,89.60,0]])

# Delta 행렬
delta = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        delta[j, i] = 2 * Bij[j, i] - Bij[i, i] - Bij[j, j]

# 2. 초기화
x = z.copy()
Psat = np.exp(A - B / (T + C))
gamma = unifgam(k, R, Q, nu, amn, n, x, T)
PHIi = np.ones(n)
K = gamma * Psat / (PHIi * P)

Aerr, Acrit = 1.0, 1e-6
alphaold = 0.5

# 3. 반복 계산
while Aerr > Acrit:
    def falpha(alpha):
        return np.sum((1 - K) * z / (1 + alpha * (K - 1)))
    
    alpha = fsolve(falpha, alphaold)[0]
    x = z / (1 + alpha * (K - 1))
    y = K * x
    
    # PHI 갱신
    gamma = unifgam(k, R, Q, nu, amn, n, x, T)
    for i in range(n):
        sumy = 0
        for j in range(n):
            for m in range(n):
                sumy += y[j] * y[m] * (2 * delta[j, i] - delta[j, m])
        PHIi[i] = np.exp((Bij[i, i] * (P - Psat[i]) + sumy * P / 2) / (Rg * T))
    
    K = gamma * Psat / (PHIi * P)
    Aerr = abs(alpha - alphaold)
    alphaold = alpha

print(f"Alpha: {alpha:.4f}")
print(f"Liquid x: {x}")
print(f"Vapor y: {y}")