import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# 1. 미분 방정식 시스템 정의
def cstrxn(X, t, Caf, D, Fi, Tf, Tj, dH, A1, E, rCp, Ui, R):
    # X[0] = h (액위), X[1] = Ca (농도), X[2] = T (온도, deg.C)
    h, Ca, T = X
    
    # 단면적(Sa) 및 열전달 면적(Sh) 계산
    Sa = (np.pi / 4) * (D**2)
    Sh = Sa + np.pi * D * h
    
    # 중간 계산값 설정[cite: 4]
    a1 = Fi / (Sa * h)
    # 온도 T를 절대온도(K)로 변환하여 속도 상수 k1 계산[cite: 4]
    k1 = A1 * np.exp(-E / (R * (T + 273.15)))
    
    b1 = dH / rCp
    c1 = (Ui * Sh) / (rCp * Sa * h)
    
    # 미분 방정식 정의 (dh/dt, dCa/dt, dT/dt)[cite: 4]
    dh = Fi / Sa - np.sqrt(10 * h / Sa)
    dCa = a1 * (Caf - Ca) - k1 * Ca
    dT = a1 * (Tf - T) + b1 * k1 * Ca - c1 * (T - Tj)
    
    return [dh, dCa, dT]

# 2. 데이터 및 파라미터 설정[cite: 4]
Caf = 10
D = 2.335
Fi = 10
Tf = 25
Tj = 25
dH = 5960
A1 = 3.49308e7
E = 11843
rCp = 500
Ui = 70
R = 1.987

# 3. 초기값 및 시간 범위 설정[cite: 4]
X0 = [1, 5, 20]  # 초기 h=1, Ca=5, T=20
t = np.linspace(0, 20, 1000)  # 0부터 20시간까지

# 4. ODE 풀기[cite: 4]
# odeint 사용 시 args를 통해 추가 파라미터를 전달합니다.
sol = odeint(cstrxn, X0, t, args=(Caf, D, Fi, Tf, Tj, dH, A1, E, rCp, Ui, R))

# 5. 결과 시각화[cite: 4]
plt.figure(figsize=(10, 8))

# h (액위) 그래프
plt.subplot(2, 2, 1)
plt.plot(t, sol[:, 0])
plt.xlabel('t(hr)')
plt.ylabel('h(m)')
plt.grid(True)

# Ca (농도) 그래프
plt.subplot(2, 2, 2)
plt.plot(t, sol[:, 1])
plt.xlabel('t(hr)')
plt.ylabel('$C_A(kmol/m^3)$')
plt.grid(True)

# T (온도) 그래프
plt.subplot(2, 2, 3)
plt.plot(t, sol[:, 2])
plt.xlabel('t(hr)')
plt.ylabel('T(deg.C)')
plt.grid(True)

plt.tight_layout()
plt.show()