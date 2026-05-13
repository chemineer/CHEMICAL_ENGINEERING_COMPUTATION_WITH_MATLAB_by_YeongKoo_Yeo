import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# 1. 파라미터 및 격자 설정
delt = 0.01
mu = 2.1e-5
C0 = 0.1
nx = 100  # 공간 격자 수
nt = 100  # 시간(또는 길이 L) 격자 수

x = np.linspace(0, 0.01, nx)
t = np.linspace(0, 1, nt)
dx = x[1] - x[0]

# 2. PDE 시스템 정의 (Method of Lines)
def pde_system(u, t):
    dudt = np.zeros_like(u)
    
    # 내부 점들에 대한 계산
    for i in range(1, nx - 1):
        # MATLAB의 c(x,t,u,DuDx) 부분: 속도 프로파일 계산
        c_val = 2 * (0.5) * (x[i]/delt - 0.5*(x[i]/delt)**2) / mu
        
        # 2차 중앙 차분 (f = DuDx이므로 df/dx = d^2u/dx^2)
        d2u_dx2 = (u[i+1] - 2*u[i] + u[i-1]) / dx**2
        
        dudt[i] = d2u_dx2 / c_val
    
    # 3. 경계 조건 (abbncon)
    # 왼쪽 (x=0): ql=1, pl=0 -> du/dx = 0 (Neumann)
    dudt[0] = (u[1] - u[0]) / dx # 단순화된 차분
    
    # 오른쪽 (x=delt): pr = ur - C0, qr=0 -> u = C0 (Dirichlet)
    # 실제 수치해석 시 경계값은 고정되므로 미분값은 0
    dudt[-1] = 0 
    
    return dudt

# 4. 초기 조건 (abitcon)
u0 = np.zeros(nx)
u0[-1] = C0  # 경계값 적용

# 5. PDE 풀기
sol = odeint(pde_system, u0, t)

# 6. 결과 시각화
fig = plt.figure(figsize=(12, 5))

# 왼쪽: 3D Surface Plot
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
X, T = np.meshgrid(x, t)
surf = ax1.plot_surface(X, T, sol, cmap='gray', edgecolor='none', antialiased=True)
ax1.set_xlabel('δ(m)')
ax1.set_ylabel('L(m)')
ax1.set_zlabel('Ca(kgmol/m³)')
ax1.set_title('Absorption Column Profile')

# 오른쪽: 2D Line Plot
ax2 = fig.add_subplot(1, 2, 2)
for k in range(0, nt, 5):  # 100개를 다 그리면 너무 겹치므로 간격을 둠
    ax2.plot(x, sol[k, :], 'k', alpha=0.5)

ax2.set_xlabel('δ(m)')
ax2.set_ylabel('Ca(kgmol/m³)')
ax2.set_title('Ca as a function of delta')

plt.tight_layout()
plt.show()