import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from pbrmf import pbrmf  # pbrmf.py 파일에서 함수를 임포트

# 1. 데이터 설정
ca0 = 1
fa0 = 1.5
k = 10
ka = 1
kb = 2
kc = 20
wf = 2

# 2. 초기 조건 및 구간 설정
wint = [0, wf]       # 독립 변수 구간 (W)
x0 = [0, 0, 0, 0]    # 초기값[cite: 1]

# 3. ODE 풀이 (MATLAB의 ode45에 해당)
# args 파라미터를 통해 추가 인자(k, fa0, ca0 등)를 전달합니다.
sol = solve_ivp(
    pbrmf, 
    wint, 
    x0, 
    args=(k, fa0, ca0, ka, kb, kc), 
    method='RK45', 
    dense_output=True
)

# 결과 추출
w = sol.t
x = sol.y  # x[0] = X1, x[1] = X2, ...

# 4. 결과 시각화[cite: 1]
plt.figure(figsize=(8, 5))
plt.plot(w, x[0], label='X_1', linestyle='-')
plt.plot(w, x[1], label='X_2', linestyle=':')
plt.plot(w, x[2], label='X_3', linestyle='-.')
plt.plot(w, x[3], label='X_4', linestyle='--')

plt.xlabel('W(kg)')
plt.ylabel('X')
plt.legend(loc='best')
plt.grid(True)
plt.show()