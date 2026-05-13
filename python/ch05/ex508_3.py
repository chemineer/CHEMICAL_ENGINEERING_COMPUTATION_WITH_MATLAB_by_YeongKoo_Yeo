import numpy as np
from scipy.optimize import fsolve
from functools import partial
from vwfun import vwfun

# 1. 파라미터 정의
T_list = np.array([40, 60, 100])
L = 1000
D = 7.981
rf = 0.00015
dz = 300
dP = -150
v0 = 10

# 결과 저장용 리스트
rhor = []
mur = []
qr = []

# 2. 온도별 루프 계산
for T in T_list:
    # 온도 T에 따른 물성치 계산 (MATLAB 식과 동일)
    rho = 62.122 + 0.0122*T - (1.54e-4)*T**2 + (2.65e-7)*T**3 - (2.24e-10)*T**4
    mu = np.exp(-11.0318 + 1057.51 / (T + 214.624))
    
    # 해당 온도와 물성치를 사용하여 vwfun 계산
    # (참고: vwfun이 내부에서 rho, mu를 매개변수로 받도록 구성되어 있다고 가정)
    func = partial(vwfun, T=T, L=L, D=D, rf=rf, dz=dz, dP=dP)
    
    # 해 찾기
    v = fsolve(func, v0)[0]
    
    # 유량 계산 (gpm)
    q = (7.481 * 60) * (np.pi * v * (D / 12)**2) / 4
    
    # 결과 저장
    rhor.append(rho)
    mur.append(mu)
    qr.append(q)

# 3. 결과 출력
print("Flow rates (qr):", qr)
print("Densities (rhor):", rhor)
print("Viscosities (mur):", mur)