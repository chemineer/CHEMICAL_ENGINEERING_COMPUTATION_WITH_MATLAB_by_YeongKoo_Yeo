import numpy as np
from phimix import phimix

# 매틀랩 코드의 변수 할당
state = 'v'
k = np.zeros((2, 2))  # 이진 상호작용 매개변수 행렬
eos = 'pr'
T = 100
P = 4.119e5
ni = np.array([0.958, 0.042]) # 조성(몰 분율)

# 임계 성질 및 편심 인자 설정
Tc = np.array([126.1, 190.6])
Pc = np.array([33.94, 46.04]) * 1e5
w = np.array([0.04, 0.011])

# 함수 호출 (Z, V, phi를 각각 반환받음)
# 매틀랩의 [Z, V, phi] = phimix(...)와 동일하게 처리
Z, V, phi = phimix(ni, P, T, Pc, Tc, w, k, state, eos)

# 결과 출력
print(f"상태 방정식({eos}) - 상태({state}) 혼합물 계산 결과:")
print(f"압축 인자 (Z) = {Z}")
print(f"몰 부피 (V) = {V}")
print(f"퓨가시티 계수 (phi) = {phi}")