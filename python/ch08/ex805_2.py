import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 기하학적 수치 설정
LA = 0.015
LB = 0.1
LC = 0.075
Lt = LA + LB + LC
kA = 0.151
kC = 0.762
qx = -15          # 열유속 (W/m^2)
T0 = 255          # x=0일 때의 초기 온도 (K)
xspan = (0, Lt)   # 위치 범위

# 2. 미분 방정식 정의 (dT/dx)
def slabmT(x, T, LA, LB, kA, kC, qx):
    # T는 배열 형태로 들어오므로 첫 번째 요소를 추출
    temp = T[0]
    
    if x <= LA:
        # A층: 일정한 열전도도 kA
        dTdx = -qx / kA
    elif x <= (LA + LB):
        # B층: 온도에 의존하는 열전도도 kB(T)
        kB = 2.5 * np.exp(-1225 / temp)
        dTdx = -qx / kB
    else:
        # C층: 일정한 열전도도 kC
        dTdx = -qx / kC
        
    return [dTdx]

# 3. ODE 풀이 (MATLAB의 ode45 대응)
# t_eval을 통해 그래프를 매끄럽게 그리기 위한 지점들을 지정합니다.
x_eval = np.linspace(0, Lt, 500)
sol = solve_ivp(
    slabmT, 
    xspan, 
    [T0], 
    args=(LA, LB, kA, kC, qx), 
    method='RK45', 
    t_eval=x_eval
)

# 4. 결과 시각화
plt.figure(figsize=(8, 5))
plt.plot(sol.t, sol.y[0])
plt.grid(True)
plt.axis([0, Lt, 250, 310])
plt.xlabel('x(m)')
plt.ylabel('T(K)')
plt.title('Temperature Profile through Multilayer Slab')
plt.show()