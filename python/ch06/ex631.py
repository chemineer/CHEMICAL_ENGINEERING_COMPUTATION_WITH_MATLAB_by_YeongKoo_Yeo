import numpy as np
from scipy.integrate import solve_ivp
from adfun import adfun  # adfun.py 임포트

# 1. 데이터 및 파라미터 설정
pf = [162, 0]           # 파라미터 벡터
Vspan = [0, 4]          # 반응기 부피 구간 (V)
X0 = [38.3, 0, 0, 1150] # 초기 조건 [FA0, FB0, FC0, T0][cite: 8]

# 2. ODE 시스템 풀이 (MATLAB의 ode45에 해당)[cite: 8]
# solve_ivp는 기본적으로 RK45 방법을 사용합니다.
sol = solve_ivp(
    adfun, 
    Vspan, 
    X0, 
    args=(pf,), 
    method='RK45', 
    dense_output=True
)

# 3. 결과 추출[cite: 8]
# sol.y는 (변수 개수, 시간 지점)의 행렬이므로 마지막 지점(-1)을 선택합니다.
FA_end = sol.y[0, -1]
FB_end = sol.y[1, -1]
FC_end = sol.y[2, -1]
T_end = sol.y[3, -1]

# 4. 결과 출력[cite: 8]
print(f'\nFA = {FA_end:g},  FB = {FB_end:g},  FC = {FC_end:g}, T  = {T_end:g}')