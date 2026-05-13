import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# 이전 단계에서 정의한 so3de 함수가 이 스크립트와 동일한 파일에 있거나 import되어야 합니다.
from so3de import so3de 

# 매개변수 설정 (MATLAB의 sa 구조체 대응)
sa = {
    'Ta': 700, 'T0': 750, 'Pt0': 202650,
    'rho0': 0.866, 'rhob': 542,
    'ya0': 0.1, 'yb0': 0.11, 'yc0': 0.79,
    'Ft0': 0.02153, 'G': 0.433,
    'epn': -0.05, 'phi': 0.45, 'mu': 3.72e-5,
    'D': 0.0453, 'Dp': 4.57e-3, 'U': 17
}

# 초기 조건 및 적분 범위
wspan = [0, 4]
z0 = [0, sa['T0'], sa['Pt0']]

# 수치 적분 (MATLAB의 ode15s에 대응하는 BDF 방법 사용)
sol = solve_ivp(
    lambda w, z: so3de(w, z, sa), 
    wspan, 
    z0, 
    method='BDF', 
    t_eval=np.linspace(wspan[0], wspan[1], 100)
)

# 결과 추출
w = sol.t
x = sol.y[0, :]
T = sol.y[1, :]
P = sol.y[2, :]

# 결과 그래프 출력
plt.figure(figsize=(10, 8))

plt.subplot(2, 2, 1)
plt.plot(w, x)
plt.grid(True)
plt.xlabel('W(kg)')
plt.ylabel('X')

plt.subplot(2, 2, 2)
plt.plot(w, T)
plt.grid(True)
plt.xlabel('W(kg)')
plt.ylabel('T(K)')

plt.subplot(2, 2, 3)
plt.plot(w, P)
plt.grid(True)
plt.xlabel('W(kg)')
plt.ylabel('P(Pa)')

plt.tight_layout()
plt.show()