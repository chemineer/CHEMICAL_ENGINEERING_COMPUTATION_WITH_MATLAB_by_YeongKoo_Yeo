import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def slab1T(x, T, A, T2, Ta, sigma):
    """
    1차원 슬래브 내부의 열전달 미분 방정식
    """
    # T는 배열 형태로 들어오므로 T[0]을 사용합니다.
    dT_dx = -sigma * (T2**4 - Ta**4) / (30 * (1 + 0.002 * T[0]) * A)
    return [dT_dx]

# 매개변수 설정
A = 1
T2_guess = 700.0  # 초기 추정치
Ta = 1273.0
sigma = 5.676e-8
x_span = [0, 0.2]
T0 = [290.0]       # 초기값 (리스트 형태)
criT = 1e-3
delT = 10.0
k = 1

# 반복법을 이용한 T2 수렴 계산
while delT > criT:
    # ODE 풀이: solve_ivp(함수, 구간, 초기값, 추가 인자)
    sol = solve_ivp(
        slab1T, 
        x_span, 
        T0, 
        args=(A, T2_guess, Ta, sigma),
        method='RK45'
    )
    
    # 마지막 지점의 온도 추출
    T_end = sol.y[0, -1]
    
    # 수렴 조건 확인 및 업데이트
    delT = abs(T2_guess - T_end)
    T2_guess = T_end
    k += 1

# 결과 출력
print(f"반복 횟수(k): {k}")
print(f"최종 온도(T_end): {T_end:.4f} K")

# 그래프 시각화
plt.figure(figsize=(8, 5))
plt.plot(sol.t, sol.y[0], label='Temperature Profile')
plt.xlabel('x(m)')
plt.ylabel('T(K)')
plt.title('One-dimensional Heat Transfer')
plt.grid(True)
plt.legend()
plt.show()