import numpy as np
from scipy.optimize import fsolve

def fd(s, F, D, z, P, Q):
    # s[0]: T0, s[1]: T1, s[2]: Tf, s[3]: T2, s[4]: T3
    # s[5]: x11, s[6]: x21, s[7]: x12, s[8]: x22, s[9]: x13, s[10]: x23
    # s[11]: V1, s[12]: V2, s[13]: V3
    
    # Antoine coefficients (1: n-butane, 2: n-pentane)
    A = np.array([6.80776, 6.85296])
    B = np.array([935.77, 1064.84])
    C = np.array([238.789, 232.012])
    
    # Enthalpy coefficients
    h1c = np.array([0.04, 29.6])
    h2c = np.array([0.025, 38.5])
    H1c = np.array([-0.04, 43.8, 8003])
    H2c = np.array([0.007, 31.7, 12004])
    
    # 평형 상수 (K) 계산 함수 (화씨 -> 섭씨 변환 포함)
    def calc_k(T):
        return 10**(A - B / ((T - 32) * 5/9 + C)) / P

    k0 = calc_k(s[0])
    k1 = calc_k(s[1])
    kf = calc_k(s[2])
    k2 = calc_k(s[3])
    k3 = calc_k(s[4])
    
    # 조성 벡터
    x1 = np.array([s[5], s[6]])
    x2 = np.array([s[7], s[8]])
    x3 = np.array([s[9], s[10]])
    x0 = k1 * x1  # 환류 액체 조성
    
    # 엔탈피 계산용 행렬 연산
    def calc_hL(T, x):
        temp_vec = np.array([T**2, T])
        coeff_mat = np.array([h1c, h2c]).T
        return np.sum((temp_vec @ coeff_mat) * x)

    def calc_hV(T, k, x):
        temp_vec = np.array([T**2, T, 1])
        coeff_mat = np.array([H1c, H2c]).T
        return np.sum((temp_vec @ coeff_mat) * (k * x))

    hL0 = calc_hL(s[0], x0)
    hL1 = calc_hL(s[1], x1)
    hLf = calc_hL(s[2], z)
    hL2 = calc_hL(s[3], x2)
    hL3 = calc_hL(s[4], x3)
    
    hV1 = calc_hV(s[1], k1, x1)
    hV2 = calc_hV(s[3], k2, x2)
    hV3 = calc_hV(s[4], k3, x3)
    
    # 유량 수지
    B_rate = F - D
    L0 = s[11] - D
    L1 = s[12] - D
    L2 = s[13] + B_rate
    L3 = B_rate
    
    # 방정식 정의 (14개)
    f = np.zeros(14)
    f[0] = np.sum(k0 * x0) - 1
    f[1] = np.sum(k1 * x1) - 1
    f[2] = np.sum(kf * z) - 1
    f[3] = np.sum(k2 * x2) - 1
    f[4] = np.sum(k3 * x3) - 1
    f[5] = -((s[11] - L0) * k1[0] + L1) * s[5] + s[12] * k2[0] * s[7]
    f[6] = -((s[11] - L0) * k1[1] + L1) * s[6] + s[12] * k2[1] * s[8]
    f[7] = -s[11] * hV1 + s[12] * hV2 - L1 * hL1 + L0 * hL0
    f[8] = L1 * s[5] - (s[12] * k2[0] + L2) * s[7] + s[13] * k3[0] * s[9] + F * z[0]
    f[9] = L1 * s[6] - (s[12] * k2[1] + L2) * s[8] + s[13] * k3[1] * s[10] + F * z[1]
    f[10] = -s[12] * hV2 + s[13] * hV3 + hLf + L1 * hL1 - L2 * hL2
    f[11] = L2 * s[7] - (s[13] * k3[0] + B_rate) * s[9]
    f[12] = L2 * s[8] - (s[13] * k3[1] + B_rate) * s[10]
    f[13] = -s[13] * hV3 + Q + L2 * hL2 - L3 * hL3
    
    return f

# 초기값 및 파라미터 설정
F = 1.0; D = 0.25; z = np.array([0.23, 0.77])
P = 760 * 120 / 14.7; Q = 1e4
s0 = [200, 145, 200, 190, 210, 0.65, 0.35, 0.43, 0.57, 0.33, 0.76, 1.1, 1, 1.1]

# 솔버 실행
s_sol = fsolve(fd, s0, args=(F, D, z, P, Q))

# 결과 출력
T0, T1, Tf, T2, T3 = s_sol[0:5]
x11, x21, x12, x22, x13, x23 = s_sol[5:11]
V1, V2, V3 = s_sol[11:14]

print(f"T0={T0:8.4f}, T1={T1:8.4f}, Tf={Tf:8.4f}, T2={T2:8.4f}, T3={T3:8.4f}")
print(f"x11={x11:6.4f}, x12={x12:6.4f}, x13={x13:6.4f}")
print(f"x21={x21:6.4f}, x22={x22:6.4f}, x23={x23:6.4f}")
print(f"V1={V1:6.4f}, V2={V2:6.4f}, V3={V3:6.4f}")