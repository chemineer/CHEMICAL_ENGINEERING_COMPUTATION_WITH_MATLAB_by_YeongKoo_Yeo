import numpy as np
from scipy.optimize import least_squares

def frx(x, t, Ca0, Ca):
    # x[0]=k’, x[1]=alpha
    k_prime = x[0]
    alpha = x[1]
    
    # MATLAB의 루프 연산을 벡터화 (7개의 잔차 반환)
    # f = x(1)*(1-x(2))*t - Ca0^(1-x(2)) + Ca(i)^(1-x(2))
    f = k_prime * (1 - alpha) * t - Ca0**(1 - alpha) + Ca**(1 - alpha)
    return f

# --- 데이터 설정 ---[cite: 6]
t = np.array([0, 50, 100, 150, 200, 250, 300], dtype=float) 
Ca = np.array([0.05, 0.038, 0.0306, 0.0256, 0.0222, 0.0195, 0.0174], dtype=float) 
Ca0 = 0.05 

# 초기 추정값: k=0.2, alpha=2[cite: 6]
initial_guess = [0.2, 2.0]

# fsolve 대신 least_squares 사용 (미지수 2개, 방정식 7개인 경우 적합)
res = least_squares(frx, initial_guess, args=(t, Ca0, Ca))

k_final, alpha_final = res.x

# 결과 출력
print(f"최적화된 속도 상수 (k'): {k_final:.6f}")
print(f"최적화된 반응 차수 (alpha): {alpha_final:.6f}")