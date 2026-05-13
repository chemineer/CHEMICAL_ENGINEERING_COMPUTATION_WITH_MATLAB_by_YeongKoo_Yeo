import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 파라미터 및 초기 설정
Cas = 0.2
R = 0.5
De = 0.1
k1a = 6.4
rspan = [0, R]
critN = 1e-6
errN = 1
Ca1 = 1e-2
Ca2 = 3e-2

# 2. 미분 방정식 정의 (rxf 함수)
def rxf(r, x, Cas, R, De, k1a):
    # x[0] = Na * r^2, x[1] = eta (유효 계수 적분항), x[2] = Ca
    if r == 0:
        Na = 0
    else:
        Na = x[0] / (r**2)
    
    # dx/dr 계산
    dx1dr = -k1a * x[2] * (r**2)
    dx2dr = 3 * x[2] * (r**2) / (Cas * (R**3))
    dx3dr = -Na / De
    
    return [dx1dr, dx2dr, dx3dr]

# 3. 이분법을 이용한 중심 농도(Ca0) 추정[cite: 13]
Cam = (Ca1 + Ca2) / 2
sol_m = None

while errN > critN:
    Cam = (Ca1 + Ca2) / 2
    
    # 각 초기값 후보에 대해 ODE 풀이 (중심 r=0에서 초기 조건 적용)[cite: 13]
    # x[0]=0 (중심에서 플럭스 0), x[1]=0 (적분 시작), x[2]=Ca_initial[cite: 13]
    sol1 = solve_ivp(rxf, rspan, [0, 0, Ca1], args=(Cas, R, De, k1a), method='RK45')
    solm = solve_ivp(rxf, rspan, [0, 0, Cam], args=(Cas, R, De, k1a), method='RK45')
    
    # 표면(r=R)에서의 농도가 Cas와 일치하는지 확인[cite: 13]
    val1 = sol1.y[2, -1] - Cas
    valm = solm.y[2, -1] - Cas
    
    if val1 * valm < 0:
        Ca2 = Cam
    else:
        Ca1 = Cam
        
    errN = abs(Ca1 - Ca2)
    sol_m = solm

# 4. 결과 계산 및 분석해 비교[cite: 13]
r = sol_m.t
Ca = sol_m.y[2, :]
effc = sol_m.y[1, -1] # 계산된 유효 계수[cite: 13]

# 분석해 계산 (Thiele modulus 이용)[cite: 13]
ephi = R * np.sqrt(k1a / De)
effa = 3 * (ephi * (1/np.tanh(ephi)) - 1) / (ephi**2)

# 결과 출력[cite: 13]
print(f'Effectiveness factor at r=R (calculated): {effc:7.5f}')
print(f'Effectiveness factor at r=R (analytic): {effa:7.5f}')

# 5. 시각화[cite: 13]
plt.figure(figsize=(8, 5))
plt.plot(r, Ca)
plt.xlabel('r(cm)')
plt.ylabel('C_A(gmol/cm^3)')
plt.title('Concentration Profile in Catalyst Pellet')
plt.grid(True)
plt.show()