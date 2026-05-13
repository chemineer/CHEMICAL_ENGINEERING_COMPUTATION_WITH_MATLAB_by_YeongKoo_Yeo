import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# 1. 데이터 및 파라미터 설정
F0 = 40; F = 40; Fj = 49.9; Ca0 = 0.55; V = 48; rho = 50; rhoj = 62.3
Cp = 0.75; Cj = 1; A = 250; U = 150
T0_val = 530; Tj0 = 530; alp = 7.08e10; lam = -3e4; E = 3e4; R = 1.9872

# 2. 온도(T)에 대한 비선형 함수 정의
def fT(T):
    # 아레니우스 항
    exp_term = alp * V * np.exp(-E / (R * T))
    
    # 에너지 수지식 (정상 상태 농도 식을 에너지 식에 대입하여 정리한 형태)
    term1 = rho * Cp * (F0 * T0_val - F * T)
    term2 = (F0 * Ca0 * V * lam * alp * np.exp(-E / (R * T))) / (F + exp_term)
    term3 = (U * A * rhoj * Cj * Fj * (T - Tj0)) / (rhoj * Cj * Fj + U * A)
    
    return term1 - term2 - term3

# 3. fsolve를 사용하여 해 구하기 (fzero 대응)
initial_guess = T0_val + 150 # 초기 추정값: 680
T_solution = fsolve(fT, initial_guess)

print(f"계산된 정상 상태 온도 (T): {T_solution[0]:.4f} R")

# 4. f(T) vs T 그래프 작성
Tv = np.arange(500, 700.1, 0.1)
Fv = fT(Tv)

plt.figure(figsize=(8, 6))
plt.plot(Tv, Fv, label='f(T)')
plt.axhline(0, color='red', linestyle='--', label='f(T) = 0') # y=0 기준선
plt.xlabel('T(R)')
plt.ylabel('f(T)')
plt.title('Nonlinear Function f(T) vs Temperature')
plt.legend()
plt.grid(True)
plt.show()