import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, bisect
from scipy.integrate import quad

# 1. 상수 및 데이터 정의
R = 83.14
w = 0.224
Tc = 304.2
Pc = 73.83
T = 288.15
Tr = T / Tc

# 2. Peng-Robinson 파라미터 계산
# alpha 및 a, b 상수 계산
kappa = 0.37464 + 1.54226 * w - 0.26992 * (w**2)
alpha = (1 + kappa * (1 - np.sqrt(Tr)))**2
a = 0.45724 * (R**2) * (Tc**2) * alpha / Pc
b = 0.0778 * R * Tc / Pc

# 3. 압력 방정식 정의 (P as a function of V)
def Ppr(V):
    return (R * T / (V - b)) - (a / (V**2 + 2 * b * V - b**2))

# 4. P vs V 그래프 그리기
V_range = np.linspace(40, 400, 300)
P_vals = Ppr(V_range)

plt.figure(figsize=(8, 5))
plt.plot(V_range, P_vals)
plt.axis([40, 400, 0, 100])
plt.grid(True)
plt.xlabel('V(cm^3/mol)')
plt.ylabel('P(bar)')
plt.show()

# 5. 등면적 법칙(Equal Area Rule)을 이용한 증기압(Pv) 계산
# 함수: P(V) - P_sat = 0, 그리고 적분 조건이 만족되는 P_sat을 찾음
# 매틀랩의 복잡한 이분법 로직을 scipy의 bisect 등으로 더 깔끔하게 처리 가능
def objective(Psat):
    # 주어진 Psat에 대해 액체 부피(Vl)와 기체 부피(Vg)를 찾음
    Vl = fsolve(lambda v: Ppr(v) - Psat, 60)[0]
    Vg = fsolve(lambda v: Ppr(v) - Psat, 300)[0]
    
    # 적분: integral(P dV) - Psat * (Vg - Vl) = 0
    integral_val, _ = quad(Ppr, Vl, Vg)
    return integral_val - Psat * (Vg - Vl)

# 이분법으로 Psat 해 찾기
Pv = bisect(objective, 40, 60)
print(f"Peng-Robinson 방정식으로 계산된 증기압 (Pv): {Pv:.4f} bar")

# 6. 확장 Antoine 방정식을 이용한 계산
Pv_antoine = 10**(47.544 - 1792.2/T - 16.559 * np.log10(T) + 0.013833 * T) / 750.0615
print(f"확장 Antoine 방정식으로 계산된 증기압 (Pv): {Pv_antoine:.4f} bar")