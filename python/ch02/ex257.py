import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 데이터 및 파라미터 설정
a = 200
k = 0.02
kg = 0.01
u = 1
Ca0 = 1
L = 1

# 초기 조건 계산 (대수 방정식 kg*(Ca - Cas) - k*Cas = 0 만족)
Cas0 = kg * Ca0 / (k + kg)
y0 = [Ca0] # Ca만 미분 방정식의 변수로 취급

# 2. 시스템 정의
# 대수 방정식을 정리하면 Cas = (kg * Ca) / (k + kg) 입니다.
# 이를 첫 번째 미분 방정식에 대입하여 단일 ODE로 만듭니다.
def hrdae_reduced(z, y):
    Ca = y[0]
    # 대수 방정식으로부터 Cas 계산
    Cas = (kg * Ca) / (k + kg)
    
    # dCa/dz 계산
    dCadz = -(kg * a / u) * (Ca - Cas)
    return [dCadz]

# 3. ODE 풀이
z_span = (0, L)
z_eval = np.linspace(0, L, 100)
sol = solve_ivp(hrdae_reduced, z_span, y0, t_eval=z_eval)

# 4. 결과 복원 (계산된 Ca를 바탕으로 Cas 재계산)
z = sol.t
Ca = sol.y[0]
Cas = (kg * Ca) / (k + kg)

# 5. 시각화
plt.figure(figsize=(8, 5))
plt.plot(z, Ca, label='$C_A$ (Fluid)')
plt.plot(z, Cas, '--', label='$C_{As}$ (Surface)')

plt.grid(True)
plt.xlabel('Location(z)')
plt.ylabel('Concentration')
plt.title('Concentration Profiles in Heterogeneous Reactor')
plt.legend()
plt.show()