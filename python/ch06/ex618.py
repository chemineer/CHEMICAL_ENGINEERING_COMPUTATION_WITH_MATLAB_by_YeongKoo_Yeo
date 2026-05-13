import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 미분 방정식 시스템 정의 (PFR 모델)
def pfrmx(t, C, k):
    # C[0] = Ch, C[1] = Cm, C[2] = Cx
    Ch, Cm, Cx = C
    
    # 반응 속도식 정의 (MATLAB 소스 기반)
    # k[0] = k1, k[1] = k2
    term1 = k[0] * (Ch**0.5) * Cm
    term2 = k[1] * Cx * (Ch**0.5)
    
    dCh = -term1 - term2
    dCm = -term1
    dCx = term1 - term2
    
    return [dCh, dCm, dCx]

# 2. 데이터 및 초기 조건 설정
k = [55.2, 30.2]          # 반응 속도 상수 k1, k2
C0 = [0.021, 0.0105, 0]   # 초기 농도 [Ch0, Cm0, Cx0]
tf = 0.5                  # 최종 체류 시간 (hr)

# 3. ODE 풀기 (solve_ivp 사용)
# t_eval을 통해 세밀한 시간 간격으로 결과를 얻어 최적점을 정밀하게 찾습니다.
t_span = (0, tf)
t_eval = np.linspace(0, tf, 1000)
sol = solve_ivp(pfrmx, t_span, C0, t_eval=t_eval, args=(k,))

# 4. 결과 데이터 추출[cite: 16]
t = sol.t
Ch = sol.y[0]
Cm = sol.y[1]
Cx = sol.y[2]

# 5. 최적 체류 시간 및 최대 농도 찾기[cite: 16]
Cxm = np.max(Cx)               # m-xylene의 최대 농도
ti = np.argmax(Cx)             # 최대 농도가 발생하는 인덱스
opmt = t[ti]                   # 최적 체류 시간

print(f"Optimum residence time = {opmt:.4f}")
print(f"Maximum concentration of m-xylene = {Cxm:.6f}")

# 6. 시각화[cite: 16]
plt.figure(figsize=(8, 6))
plt.plot(t, Ch, '--', label='$C_H$ (Hydrogen)')
plt.plot(t, Cm, ':', label='$C_M$ (Mesitylene)')
plt.plot(t, Cx, '-', label='$C_X$ (m-Xylene)')

# 최적점 표시
plt.plot(opmt, Cxm, 'ro', label='Optimum Point')

plt.xlabel('$\\tau$(hr)')
plt.ylabel('$C(lbmol/ft^3)$')
plt.title('Concentration Profiles in PFR')
plt.legend()
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()