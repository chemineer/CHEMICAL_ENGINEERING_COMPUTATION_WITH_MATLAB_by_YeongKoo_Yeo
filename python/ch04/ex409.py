import numpy as np
from khmix import khmix

# 1. 초기값 설정
state = 'L'
nx = np.array([0.419, 0.3783, 0.2027]) # 액체 조성
P = 6.8947e5
T = 158
eos = 'rk'

# 2. 물성치 및 파라미터 설정
Pc = np.array([45.99, 48.72, 42.48]) * 1e5
Tc = np.array([190.6, 305.3, 369.8])
w = np.array([0.012, 0.1, 0.152])
k = np.zeros((3, 3))

# 열용량 계수 행렬 (Afi)
Afi = np.array([
    [8.245223, 0.3806333e-2, 0.8864745e-5, -0.7461153e-8, 0.182296e-11],
    [11.51606, 0.140309e-1, 0.854034e-5, -0.1106078e-7, 0.31622e-11],
    [15.58683, 0.2504953e-1, 0.1404258e-4, -0.352626e-7, 0.1864467e-10]
])

# 3. 상태 결정 (액체일 경우 nx 사용)
state_upper = state.upper()
if state_upper == 'L':
    x = nx
else:
    # ny가 정의되어 있지 않으므로 기체일 경우의 처리가 필요합니다.
    x = None 

# 4. 혼합물 물성 계산 함수 호출
Z, H = khmix(x, P, T, state_upper, eos, Pc, Tc, w, k, Afi)

# 5. 결과 출력
print(f"Equation of state: {eos.upper()}, State: {state_upper}")
print(f"Compressibility factor = {Z}, Enthalpy = {H} (J/mol)")