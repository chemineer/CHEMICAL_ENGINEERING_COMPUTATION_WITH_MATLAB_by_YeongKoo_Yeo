import numpy as np
from scipy.optimize import fsolve
from pipnet import pipnet

# 데이터 설정
rho = 881
mu = 5e-4
D = np.array([0.508, 0.406, 0.406, 0.610, 0.508, 0.406, 0.305, 0.406, 0.305, 0.406, 0.305, 0.305])
L = np.array([915, 1220, 915, 1220, 915, 1220, 915, 1220, 1220, 915, 915, 1220])
rf = 5e-6 * np.ones(len(D))
m0 = 20 * np.ones(len(D))

# fsolve 호출
# args를 통해 pipnet 함수에 필요한 추가 인자를 전달합니다.
m = fsolve(pipnet, m0, args=(rho, mu, D, L, rf))

print("계산된 질량 유량 (m):", m)