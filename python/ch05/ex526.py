import numpy as np
# dpcatbed.py 모듈로부터 dpcatbed 함수를 임포트합니다.
from dpcatbed import dpcatbed

# --- 입력 데이터 설정 (source: 4) ---
Bd = 0.134     # Bed diameter
Bl = 60.0      # Bed length
Pd = 0.25      # Particle diameter[cite: 4]
Pl = 0.25      # Particle length[cite: 4]
ep = 0.45      # Void fraction (epsilon)[cite: 4]
W = 104.4      # Mass flow rate[cite: 4]
mu = 0.0278    # Viscosity[cite: 4]
rho = 0.413    # Density[cite: 4]

# --- 충전층 압력 강하 계산 ---
# MATLAB: dPt = dpcatbed(W,Bd,Bl,Pd,Pl,ep,mu,rho)[cite: 4]
dPt = dpcatbed(W, Bd, Bl, Pd, Pl, ep, mu, rho)

# --- 결과 출력 ---
print(f"Total Pressure Drop (dPt): {dPt:.6f}")