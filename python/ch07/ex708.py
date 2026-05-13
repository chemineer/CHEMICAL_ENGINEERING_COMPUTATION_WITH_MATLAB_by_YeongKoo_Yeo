import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 파라미터 설정
Dab = 1.5e-9
delt = 3e-4
delx = delt / 10
vm = 0.6
cas = 0.03
kp = 0
zf = 1
zspan = (0, zf)

# 2. 초기 조건 설정 (x(1) ~ x(9) 모두 0)
c0 = np.zeros(9)

# 3. 미분 방정식 정의 (ODE function)
def bsf(z, x, Dab, kp, delx, delt, vm, cas):
    # x[0]=Ca2, x[1]=Ca3, ..., x[8]=Ca10
    ca1 = cas
    
    # 하류 경계 조건 (Boundary Condition)
    if (4 * x[8] < x[7]):
        ca11 = 0
    else:
        ca11 = (4 * x[8] - x[7]) / 3
        
    dxdz = np.zeros(9)
    # 각 노드에 대한 변화율 계산
    # MATLAB의 (1*delx/delt) 등 인덱스 기반 속도 분포 반영
    for i in range(9):
        # i=0일 때 이전 값은 ca1, i=8일 때 다음 값은 ca11
        prev_x = ca1 if i == 0 else x[i-1]
        next_x = ca11 if i == 8 else x[i+1]
        
        # 속도 항 (Velocity profile) 계산
        v_z = vm * (1 - ((i + 1) * delx / delt)**2)
        
        # 확산 및 반응 항 계산
        dxdz[i] = (Dab * (next_x - 2 * x[i] + prev_x) / (delx**2) - kp * x[i]) / v_z
        
    return dxdz

# 4. ODE 풀이
sol = solve_ivp(
    bsf, 
    zspan, 
    c0, 
    args=(Dab, kp, delx, delt, vm, cas), 
    method='RK45',
    t_eval=np.linspace(0, zf, 500) # 그래프 출력을 위한 지점 설정
)

z = sol.t
x = sol.y.T # (z_points, species_points)

# 5. 경계값(Ca1, Ca11) 복원 및 전체 데이터 구성
nz = len(z)
xc1 = cas * np.ones(nz)
xc11 = np.zeros(nz)

for k in range(nz):
    if (4 * x[k, 8] < x[k, 7]):
        xc11[k] = 0
    else:
        xc11[k] = (4 * x[k, 8] - x[k, 7]) / 3

# Ca1 ~ Ca11 합치기
c = np.column_stack((xc1, x, xc11))

# 6. 결과 시각화
plt.figure(figsize=(10, 6))
# 파이썬 인덱스 2, 4, 6, 8은 MATLAB의 c(:, 3, 5, 7, 9)와 대응
plt.plot(z, c[:, 2], label='$C_{A3}$')
plt.plot(z, c[:, 4], ':', label='$C_{A5}$')
plt.plot(z, c[:, 6], '-.', label='$C_{A7}$')
plt.plot(z, c[:, 8], '--', label='$C_{A9}$')

plt.xlabel('z(m)')
plt.ylabel('C(kgmol/m$^3$)')
plt.legend(loc='best')
plt.grid(True, alpha=0.3)
plt.title('Concentration Profiles along z')
plt.show()