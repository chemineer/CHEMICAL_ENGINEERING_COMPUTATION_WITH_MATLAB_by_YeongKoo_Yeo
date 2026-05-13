import numpy as np
import matplotlib.pyplot as plt

# 1. 데이터 및 파라미터 초기화
A1, A2 = 1.1, 172.2
E1, E2 = 2.09e4, 4.18e4
R = 8.314
dH1, dH2 = 4.18e4, 8.36e4
rho, Cp = 1000, 1
Tc = 25
Tsmax, Tsmin = 150, 70
Uj = 1.16
Ucmax, Ucmin = 4.42, 1.39
AcV, AjV = 17, 30
Kc = 0.1
tauI = 360
us = 1
dt = 0.1
tmax = 4000

# 파라미터 계산
gam1 = dH1 / (rho * Cp)
gam2 = dH2 / (rho * Cp)
a1 = (Uj * Tsmin * AjV + Ucmax * Tc * AcV) / (rho * Cp)
a2 = -(Uj * AjV + Ucmax * AcV) / (rho * Cp)
b1 = (Uj * AjV * (Tsmax - Tsmin) - (Ucmax - Ucmin) * AcV * Tc) / (rho * Cp)
b2 = (Ucmax - Ucmin) * AcV / (rho * Cp)

# 2. 초기화 및 배열 생성
t = np.arange(0, tmax + dt * 1.5, dt)
n = len(t)
T = np.zeros(n); Ca = np.zeros(n); Cb = np.zeros(n)
u = np.zeros(n); er = np.zeros(n); erc = np.zeros(n)
Ts = np.zeros(n); Uc = np.zeros(n); Fc = np.zeros(n)

T[0], Ca[0], Cb[0] = 25, 1, 0
Td = 54 + 71 * np.exp(-0.0025 * t)  # 목표 온도 궤적

er[0] = 100
erc[0] = 0
u[0] = us
Ts[0] = (Tsmax - Tsmin) * u[0] + Tsmin
Uc[0] = (Ucmin - Ucmax) * u[0] + Ucmax
Fc[0] = (4550 * (1 / Uc[0] - 1 / 10.8))**(-1.25)

# 3. 미분 방정식 정의 (ODE 함수)
def ode_system(T_val, Ca_val, Cb_val, u_val):
    temp_k = 273.15 + T_val
    rate1 = A1 * np.exp(-E1 / (R * temp_k)) * Ca_val**2
    rate2 = A2 * np.exp(-E2 / (R * temp_k)) * Cb_val
    
    dTdt = gam1 * rate1 + gam2 * rate2 + (a1 + a2 * T_val) + (b1 + b2 * T_val) * u_val
    dCadt = -rate1
    dCbdt = rate1 - rate2
    return dTdt, dCadt, dCbdt

# 4. 시뮬레이션 루프 (RK4 + PI Control)
for k in range(n - 1):
    # Runge-Kutta 4th Order
    k1, k11, k12 = ode_system(T[k], Ca[k], Cb[k], u[k])
    
    k2, k21, k22 = ode_system(T[k] + k1*dt/2, Ca[k] + k11*dt/2, Cb[k] + k12*dt/2, u[k])
    
    k3, k31, k32 = ode_system(T[k] + k2*dt/2, Ca[k] + k21*dt/2, Cb[k] + k22*dt/2, u[k])
    
    k4, k41, k42 = ode_system(T[k] + k3*dt, Ca[k] + k31*dt, Cb[k] + k32*dt, u[k])
    
    T[k+1] = T[k] + dt * (k1/6 + k2/3 + k3/3 + k4/6)
    Ca[k+1] = Ca[k] + dt * (k11/6 + k21/3 + k31/3 + k41/6)
    Cb[k+1] = Cb[k] + dt * (k12/6 + k22/3 + k32/3 + k42/6)
    
    # PI 제어 및 Anti-windup (Saturation)
    erc[k+1] = erc[k] + er[k] * dt
    er[k+1] = Td[k+1] - T[k+1]
    u[k+1] = us + Kc * (er[k+1] + erc[k+1] / tauI)
    
    # u값 제한 (0 <= u <= 1)
    u[k+1] = np.clip(u[k+1], 0, 1)
    
    # 파라미터 업데이트
    Ts[k+1] = (Tsmax - Tsmin) * u[k+1] + Tsmin
    Uc[k+1] = (Ucmin - Ucmax) * u[k+1] + Ucmax
    Fc[k+1] = (4550 * (1 / Uc[k+1] - 1 / 10.8))**(-1.25)

# 5. 결과 시각화
plt.figure(figsize=(12, 5))

# Subplot 1: 농도 변화
plt.subplot(1, 2, 1)
plt.plot(t, Ca, '--', label='$C_A(t)$')
plt.plot(t, Cb, label='$C_B(t)$')
plt.xlabel('t(sec)')
plt.ylabel('C(kmol/m^3)')
plt.legend()
plt.axis([0, 4000, 0, 1])
plt.grid(True)

# Subplot 2: 온도 변화
plt.subplot(1, 2, 2)
plt.plot(t, Td, '--', label='Desired temp.')
plt.plot(t, T, label='Reactor temp.')
plt.plot(t, Ts, '.-', label='Steam temp.', markevery=200) # 점이 너무 많아 간격 조절
plt.xlabel('t(sec)')
plt.ylabel('Temp(deg.C)')
plt.legend()
plt.axis([0, 4000, 20, 160])
plt.grid(True)

plt.tight_layout()
plt.show()