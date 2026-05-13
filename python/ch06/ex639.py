import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# 1. 데이터 구조 설정 (클래스 또는 딕셔너리)
class MBData:
    def __init__(self):
        self.Y = 0.4
        self.D = 0.202
        self.alp1 = 2.2
        self.alp2 = 0.2
        self.Pf = 50
        self.mum = 0.48
        self.K1 = 0.04545
        self.Km = 1.2
        self.Sf = 20

mbdat = MBData()

# 2. 미분 방정식 정의 (mbrxn 함수)
def mbrxn(t, z, mbdat):
    # z[0]=x, z[1]=S, z[2]=P
    x, S, P = z
    
    # 비성장 속도 mu 계산 (Haldane 식에 생성물 저해 효과 포함)
    mu = mbdat.mum * S * (1 - P / mbdat.Pf) / (mbdat.Km + S + mbdat.K1 * (S**2))
    
    # 미분값 계산[cite: 18]
    dxdt = (mu - mbdat.D) * x
    dSdt = mbdat.D * (mbdat.Sf - S) - mu * x / mbdat.Y
    dPdt = -mbdat.D * P + (mbdat.alp1 * mu + mbdat.alp2) * x
    
    return [dxdt, dSdt, dPdt]

# 3. 정상 상태 방정식 정의 (mbsrxn 함수)[cite: 18]
def mbsrxn(z, mbdat):
    # 시간 t가 없는 형태의 미분 방정식 반환 (fsolve용)[cite: 18]
    return mbrxn(0, z, mbdat)

# 4. ODE 풀이 (동적 거동 분석)[cite: 18]
x0, S0, P0 = 1, 50, 0
z0 = [x0, S0, P0]
tint = [0, 120]

sol = solve_ivp(mbrxn, tint, z0, args=(mbdat,), method='RK45', t_eval=np.linspace(0, 120, 500))

t = sol.t
x_val, S_val, P_val = sol.y

# 5. 정상 상태 값 계산 (fsolve)[cite: 18]
zs0 = [5, 5, 10] # 초기 추정치[cite: 18]
zs = fsolve(mbsrxn, zs0, args=(mbdat,))

print(f'At steady-state, x = {zs[0]:g}, S = {zs[1]:g}, P = {zs[2]:g}')

# 6. 결과 시각화[cite: 18]
plt.figure(figsize=(12, 5))

# Subplot 1: x(세포)와 P(생성물) 농도[cite: 18]
plt.subplot(1, 2, 1)
plt.plot(t, x_val, label='x')
plt.plot(t, P_val, '--', label='P')
plt.grid(True)
plt.xlabel('t(h)')
plt.ylabel('x and P (g/l)')
plt.legend(loc='best')

# Subplot 2: S(기질) 농도[cite: 18]
plt.subplot(1, 2, 2)
plt.plot(t, S_val)
plt.grid(True)
plt.xlabel('t(h)')
plt.ylabel('S (g/l)')

plt.tight_layout()
plt.show()