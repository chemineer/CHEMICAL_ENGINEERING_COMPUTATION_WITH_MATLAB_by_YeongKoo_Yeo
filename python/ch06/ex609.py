import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# 1. 데이터 및 파라미터 설정
Ca0 = 0.4
k0 = 460
E = 1380
tau = 0.18
Tc = 298.15
kappa = 78
Cp = 32

# 2. 미분 방정식 정의 (dy/dt)
def dy_dt(y, t):
    Ca = y[0]
    T = y[1]
    
    # 반응 속도 상수 k
    k = k0 * np.exp(-E / T)
    
    # 농도 변화율 dCa/dt
    dCa = (Ca0 - Ca) / tau - (k * Ca)
    
    # 온도 변화율 dT/dt
    # MATLAB 소스의 이중 음수(--) 및 복잡한 식을 파이썬에 맞게 정리
    dT = ((-(-151080 + 2 * (T - 298.15)) / Cp) * (k * Ca / Ca0) - 
          (1 + kappa) * (T - Tc) / tau)
    
    return [dCa, dT]

# 3. 초기 조건 및 시간 범위 설정[cite: 2]
y0 = [0.1, 300]  # 초기 농도 0.1, 초기 온도 300K
t = np.linspace(0, 1, 500)  # 0부터 1시간까지 500개 지점

# 4. ODE 풀기[cite: 2]
sol = odeint(dy_dt, y0, t)

# 5. 결과 가공 (K -> deg.C 변환)[cite: 2]
Ca_result = sol[:, 0]
T_degC = sol[:, 1] - 273.15

# 6. 시각화[cite: 2]
plt.figure(figsize=(12, 5))

# 농도 그래프
plt.subplot(1, 2, 1)
plt.plot(t, Ca_result, 'b-')
plt.xlabel('t(hr)')
plt.ylabel('$C_A(mol/cm^3)$')
plt.title('Concentration of A')
plt.grid(True)

# 온도 그래프
plt.subplot(1, 2, 2)
plt.plot(t, T_degC, 'r-')
plt.xlabel('t(hr)')
plt.ylabel('T(deg.C)')
plt.title('Temperature')
plt.grid(True)

plt.tight_layout()
plt.show()