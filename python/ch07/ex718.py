import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 구조 정의 (Dictionary 사용)
cmdat = {
    'F': 100, 'z': 0.5, 'R': 128.01, 'V': 178.01, 'mr': 10,
    'mf': 10, 'ms': 10, 'md': 100, 'mb': 100, 'alpa': 2
}

# 2. ODE 시스템 정의 (cmeqn)
def cmeqn(t, x, cmdat):
    # x(0)=xb, x(1)=xs, x(2)=xf, x(3)=xr, x(4)=xd
    F, z, R, V, alpa = cmdat['F'], cmdat['z'], cmdat['R'], cmdat['V'], cmdat['alpa']
    mr, mf, ms, md, mb = cmdat['mr'], cmdat['mf'], cmdat['ms'], cmdat['md'], cmdat['mb']
    
    # 평형 관계식 (Equilibrium)
    def get_y(xi):
        return alpa * xi / (1 + (alpa - 1) * xi)
    
    yb = get_y(x[0])
    ys = get_y(x[1])
    yf = get_y(x[2])
    yr = get_y(x[3])
    
    # 유량 정의
    Lr = R
    Ls = R + F
    B = Ls - V
    
    # 미분 방정식 (dx/dt)
    dxb = (Ls * x[1] - B * x[0] - V * yb) / mb
    dxs = (Ls * (x[2] - x[1]) + V * (yb - ys)) / ms
    dxf = (Lr * (x[3] - x[2]) + F * (z - x[2]) + V * (ys - yf)) / mf
    dxr = (Lr * (x[4] - x[3]) + V * (yf - yr)) / mr
    dxd = V * (yr - x[4]) / md
    
    return [dxb, dxs, dxf, dxr, dxd]

# 3. 초기값 및 시간 범위 설정
x0 = [0, 0, cmdat['z'], 0, 0]  # xb0, xs0, xf0, xr0, xd0
t_span = (0, 12)
t_eval = np.linspace(0, 12, 300) # 매끄러운 그래프를 위한 지점

# 4. ODE 풀기 (solve_ivp)
sol = solve_ivp(cmeqn, t_span, x0, args=(cmdat,), t_eval=t_eval, method='RK45')

# 5. 결과 시각화
plt.figure(figsize=(10, 6))
plt.plot(sol.t, sol.y[0], label='x_B')
plt.plot(sol.t, sol.y[1], '--', label='x_S')
plt.plot(sol.t, sol.y[2], ':', label='x_f')
plt.plot(sol.t, sol.y[3], '-.', label='x_R')
plt.plot(sol.t, sol.y[4], '.', label='x_D', markersize=4)

plt.grid(True)
plt.xlabel('t(min)')
plt.ylabel('Composition')
plt.legend(loc='best')
plt.title('Dynamic Composition of Distillation Column')
plt.show()