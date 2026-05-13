# main.py
import numpy as np
from scipy.optimize import fsolve
from cstrmult import cstrmult  # 제작한 모듈 import

# 1. 데이터 및 초기 파라미터 설정
ka = 10
kc = 15
V = 2500
v0 = 100
Ca0 = 2
Cb0 = 2
C0 = [Ca0, Cb0, 0, 0]  # 초기 추정값[cite: 8]

# 2. fsolve를 사용하여 비선형 방정식 해 구하기[cite: 8]
# args를 통해 추가 파라미터를 전달합니다.
C_sol = fsolve(cstrmult, C0, args=(ka, kc, Ca0, Cb0, v0, V))

Ca_f, Cb_f, Cc_f, Cd_f = C_sol

# 3. 선택도(Selectivity) 계산[cite: 8]
# Cd가 매우 작을 경우 0으로 처리하는 로직 포함[cite: 8]
if Cd_f <= 1e-4:
    Scd_f = 0
else:
    Scd_f = Cc_f / Cd_f

# 4. 결과 출력[cite: 8]
print('The final concentration of each species:')
print(f'Caf = {Ca_f:.6f}, Cbf = {Cb_f:.6f}, Ccf = {Cc_f:.6f}, Cdf = {Cd_f:.6f}')
print(f'The final selectivity: Scdf = {Scd_f:.6f}')