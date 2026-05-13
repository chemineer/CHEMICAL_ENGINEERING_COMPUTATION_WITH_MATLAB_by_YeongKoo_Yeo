import numpy as np
from scipy.optimize import fsolve

def debutanizer(x, zf, q, relvol, d_base):
    # x(0): theta1, x(1): theta2, x(2): d(4), x(3): D, x(4): Rmin
    # 파이썬은 인덱스가 0부터 시작함에 유의하세요.
    
    d = np.array(d_base, dtype=float)
    d[3] = x[2]  # MATLAB의 d(4) = x(3)
    
    theta1 = x[0]
    theta2 = x[1]
    d_comp_4 = x[2]
    D = x[3]
    Rmin = x[4]
    
    # 방정식 정의
    f1 = np.sum(d) - D
    f2 = np.sum(relvol * zf / (relvol - theta1)) - (1 - q)
    f3 = np.sum(relvol * zf / (relvol - theta2)) - (1 - q)
    f4 = np.sum(relvol * d / (relvol - theta1)) - D * (1 + Rmin)
    f5 = np.sum(relvol * d / (relvol - theta2)) - D * (1 + Rmin)
    
    return [f1, f2, f3, f4, f5]

# 초기 데이터 설정
zf = np.array([0.0137, 0.5113, 0.0411, 0.0171, 0.0262, 0.0446, 0.3106, 0.0354])
relvol = np.array([2.43, 1.93, 1.00, 0.765, 0.362, 0.164, 0.0720, 0.0362])
d_base = [12, 442, 13, 0, 0, 0, 0, 0]
q = 0.87
x0 = [1.5, 0.8, 3, 500, 1.5]

# 비선형 방정식 풀이
sol = fsolve(debutanizer, x0, args=(zf, q, relvol, d_base))

# 결과 출력
print("--- 계산 결과 ---")
print(f"theta 1: {sol[0]:.4f}")
print(f"theta 2: {sol[1]:.4f}")
print(f"d(4)   : {sol[2]:.4f}")
print(f"D      : {sol[3]:.4f}")
print(f"Rmin   : {sol[4]:.4f}")