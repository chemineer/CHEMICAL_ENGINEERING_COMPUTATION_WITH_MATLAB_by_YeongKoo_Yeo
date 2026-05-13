import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# 1. 데이터 및 파라미터 설정
Caf = 10; D = 2.335; Fi = 10; Tf = 25; Tj = 25; dH = 5960; A1 = 3.49308e7
E = 11843; rCp = 500; Ui = 70; R = 1.987

# 2. 고정된 액위(hs) 및 관련 계수 계산
Sa = (np.pi / 4) * D**2
hs = Fi**2 / (10 * Sa)
Sh = Sa + np.pi * D * hs
a1 = Fi / Sa / hs
b1 = dH / rCp
c1 = Ui * Sh / (rCp * Sa * hs)

# 3. 정상 상태 방정식 정의
def sscstr(X):
    Cas, Ts = X  # X[0] = 농도, X[1] = 온도(deg.C)
    k1 = A1 * np.exp(-E / (R * (Ts + 273.15)))
    
    # f1: 농도 수지, f2: 에너지 수지
    f1 = a1 * (Caf - Cas) - k1 * Cas
    f2 = a1 * (Tf - Ts) + b1 * k1 * Cas - c1 * (Ts - Tj)
    return [f1, f2]

# 4. 세 가지 서로 다른 초기값으로 해 구하기
guesses = [[8, 20], [5, 70], [2, 120]]
solutions = []

print("CSTR Steady-State Analysis:")
for x0 in guesses:
    sol = fsolve(sscstr, x0)
    solutions.append(sol)
    print(f"Initial guess: Cas = {x0[0]:g}, Ts = {x0[1]:g}", end="")
    print(f"  -> Steady-state: Cas = {sol[0]:.4f}, Ts = {sol[1]:.2f}")

# 5. 열 및 온도 프로필 생성 (그래프용 데이터)[cite: 5]
Ts_range = np.arange(20, 120.1, 0.1)
k1s = A1 * np.exp(-E / (R * (Ts_range + 273.15)))
Cas_range = a1 * Caf / (k1s + a1)

# 제거되는 열(Qr)과 발생되는 열(Qg)[cite: 5]
Qr = Ui * Sh * (Ts_range - Tj) + Fi * rCp * (Ts_range - Tf)
Qg = dH * Sa * hs * Cas_range * k1s

# 자켓 온도(Tjs) 프로필 계산[cite: 5]
Tjs = Ts_range + (rCp * Fi * (Ts_range - Tf) - dH * Sa * hs * Cas_range * k1s) / (Ui * Sh)

# 6. 결과 시각화[cite: 5]
plt.figure(figsize=(12, 5))

# 왼쪽 그래프: 열 프로필 (S-곡선 확인)[cite: 5]
plt.subplot(1, 2, 1)
plt.plot(Ts_range, Qr, ':', label='Qr (Heat Removal)')
plt.plot(Ts_range, Qg, '-', label='Qg (Heat Generation)')
plt.xlabel('Reactor temp.(deg.C)')
plt.ylabel('Q(kcal/h)')
plt.title('Heat Generation vs removal')
plt.legend()
plt.grid(True)

# 오른쪽 그래프: Tj에 따른 Reactor Temperature[cite: 5]
plt.subplot(1, 2, 2)
plt.plot(Tjs, Ts_range)
plt.axis([0, 60, 25, 110])
plt.xlabel('$T_j$(deg.C)')
plt.ylabel('T(deg.C)')
plt.title('Hysteresis Loop (T vs Tj)')
plt.grid(True)

plt.tight_layout()
plt.show()