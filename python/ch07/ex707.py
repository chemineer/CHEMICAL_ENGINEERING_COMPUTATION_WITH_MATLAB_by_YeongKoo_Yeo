import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 파라미터 설정
Dab = 1e-9
delx = 5e-4
ca10 = 1e-3
ca90 = 2e-3
ca0 = 6e-3
K = 1.5
tf = 20000
tspan = (0, tf)

# 2. 초기 조건 설정 (k=1 to 9 -> x0[0] to x0[8])
# MATLAB의 x0(2) ~ x0(8)이 파이썬의 초기값 리스트가 됩니다.
x0_full = [ca10 + (ca90 - ca10) * (k) / 8 for k in range(9)]
c0 = x0_full[1:8]  # x0(2)부터 x0(8)까지 추출

# 3. 미분 방정식 정의 (ODE function)
def slf(t, x, Dab, delx, ca0, ca10, ca90, K):
    # x[0]=Ca2, x[1]=Ca3, ..., x[6]=Ca8
    if t == 0:
        ca1 = ca10
        ca9 = ca90
    else:
        ca1 = ca0 / K
        ca9 = (4 * x[6] - x[5]) / 3
    
    # 각 노드별 변화량 계산
    dx = np.zeros(7)
    dx[0] = Dab * (x[1] - 2 * x[0] + ca1) / (delx**2)
    dx[1] = Dab * (x[2] - 2 * x[1] + x[0]) / (delx**2)
    dx[2] = Dab * (x[3] - 2 * x[2] + x[1]) / (delx**2)
    dx[3] = Dab * (x[4] - 2 * x[3] + x[2]) / (delx**2)
    dx[4] = Dab * (x[5] - 2 * x[4] + x[3]) / (delx**2)
    dx[5] = Dab * (x[6] - 2 * x[5] + x[4]) / (delx**2)
    dx[6] = Dab * (ca9 - 2 * x[6] + x[5]) / (delx**2)
    
    return dx

# 4. ODE 풀이 (ode45 대신 solve_ivp 사용)
sol = solve_ivp(
    slf, 
    tspan, 
    c0, 
    args=(Dab, delx, ca0, ca10, ca90, K), 
    method='RK45', 
    t_eval=np.linspace(0, tf, 500) # 그래프를 매끄럽게 하기 위해 시간축 설정
)

t = sol.t
x = sol.y.T # (시간, 변수) 형태로 전치

# 5. 경계값 계산 (Ca1 및 Ca9 계산)
nt = len(t)
# Ca1 계산
xc1 = np.full(nt, ca0 / K)
xc1[0] = ca10

# Ca9 계산
xc9 = (4 * x[:, 6] - x[:, 5]) / 3
xc9[0] = ca90

# 전체 데이터 합치기 (Ca1 ~ Ca9)
c = np.column_stack((xc1, x, xc9))

# 6. 결과 시각화
plt.figure(figsize=(10, 6))
plt.plot(t, c[:, 2], label='$C_{A3}$')       # c[:, 2]는 MATLAB의 c(:, 3)
plt.plot(t, c[:, 4], ':', label='$C_{A5}$')    # c[:, 4]는 MATLAB의 c(:, 5)
plt.plot(t, c[:, 6], '-.', label='$C_{A7}$')   # c[:, 6]는 MATLAB의 c(:, 7)
plt.plot(t, c[:, 8], '--', label='$C_{A9}$')   # c[:, 8]는 MATLAB의 c(:, 9)

plt.xlabel('t(sec)')
plt.ylabel('C(kgmol/m$^3$)')
plt.legend(loc='best')
plt.grid(True)
plt.title('Diffusion Profile Over Time')
plt.show()