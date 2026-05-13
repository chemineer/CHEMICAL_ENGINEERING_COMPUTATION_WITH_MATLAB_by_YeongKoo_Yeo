import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 1. 파라미터 및 데이터 설정
class Parameter:
    def __init__(self):
        self.n = 50          # 격자 수
        self.L = 2           # 반응기 길이 (m)
        self.v = 0.4         # 유속 (m/min)
        self.Cf = 1          # 유입 농도 (mol/m^3)
        self.Tf = 450        # 유입 온도 (K)
        self.E = 6e4         # 활성화 에너지 (J/mol)
        self.dEr = -1e4      # 흡착 에너지 변화량
        self.dH = -1e5       # 반응열 (J/mol)
        self.rCp = 800       # 밀도 * 비열 (J/m^3·K)
        self.T1 = 450        # 기준 온도 (K)
        self.k1 = 0.2        # 기준 속도 상수
        self.Kr1 = 1         # 기준 흡착 상수
        self.R = 8.314       # 기체 상수

pf = Parameter()

# 2. 미분 방정식 시스템 정의 (Method of Lines)
def dfr(t, Z, pf):
    # Z는 [C1, T1, C2, T2, ..., Cn, Tn] 형태의 1차원 배열
    Z_reshaped = Z.reshape((pf.n, 2))
    C = Z_reshaped[:, 0]
    T = Z_reshaped[:, 1]
    
    dC = np.zeros(pf.n)
    dT = np.zeros(pf.n)
    h = pf.L / pf.n
    
    for i in range(pf.n):
        # 아레니우스 식을 이용한 온도 의존적 상수 계산
        k = pf.k1 * np.exp(-(pf.E / pf.R) * (1/T[i] - 1/pf.T1))
        Kr = pf.Kr1 * np.exp(-(pf.dEr / pf.R) * (1/T[i] - 1/pf.T1))
        
        # 반응 속도식 (Langmuir-Hinshelwood 형태 가정)
        rx = k * C[i] / np.sqrt(1 + Kr * C[i]**2)
        
        # 상류(Upwind) 차분법 적용[cite: 18]
        if i == 0:
            s = (pf.v / h) * (C[i] - pf.Cf)
            d = (pf.v / h) * (T[i] - pf.Tf)
        else:
            s = (pf.v / h) * (C[i] - C[i-1])
            d = (pf.v / h) * (T[i] - T[i-1])
        
        # 물질 수지 및 에너지 수지[cite: 18]
        dC[i] = -s - rx
        dT[i] = -d + (-pf.dH) * rx / pf.rCp
    
    # 다시 1차원 배열로 반환[cite: 18]
    return np.column_stack((dC, dT)).flatten()

# 3. 초기 조건 및 시간 설정[cite: 18]
h = pf.L / pf.n
w = np.arange(1, pf.n + 1) * h  # 공간 그리드 (x축)
# 초기 상태: 모든 격자점이 유입 농도와 온도로 채워져 있다고 가정[cite: 18]
Z0 = np.tile([pf.Cf, pf.Tf], pf.n)
t_span = (0, 10)
t_eval = np.arange(0, 10.1, 0.1)

# 4. ODE 풀기 (강성 시스템을 위해 BDF 솔버 사용)[cite: 18]
sol = solve_ivp(dfr, t_span, Z0, method='BDF', t_eval=t_eval, args=(pf,))

# 5. 결과 재구성[cite: 18]
# sol.y의 형상은 (2*n, t_steps)이므로 전치하여 사용
Z_res = sol.y.T
C = Z_res[:, 0::2] # 0부터 2씩 건너뛰며 추출 (농도)
T = Z_res[:, 1::2] # 1부터 2씩 건너뛰며 추출 (온도)

# 6. 시각화 (3D Mesh Plot)[cite: 18]
fig = plt.figure(figsize=(14, 6))
W, T_mesh = np.meshgrid(w, sol.t)

# 농도 분포
ax1 = fig.add_subplot(1, 2, 1, projection='3d')
surf1 = ax1.plot_surface(W, T_mesh, C, cmap='gray', edgecolor='none', alpha=0.8)
ax1.set_xlabel('x (m)')
ax1.set_ylabel('t (min)')
ax1.set_zlabel('C (mol/m^3)')
ax1.set_title('Concentration Profile')

# 온도 분포
ax2 = fig.add_subplot(1, 2, 2, projection='3d')
surf2 = ax2.plot_surface(W, T_mesh, T, cmap='gray', edgecolor='none', alpha=0.8)
ax2.set_xlabel('x (m)')
ax2.set_ylabel('t (min)')
ax2.set_zlabel('T (K)')
ax2.set_title('Temperature Profile')

plt.tight_layout()
plt.show()