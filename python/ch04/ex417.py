import numpy as np
from scipy.optimize import fsolve
from binflash import binflash

# 1. 입력 데이터 정의
v = 0.65       # 기화율
T = 60.0       # 온도 (C)
z = 0.6        # 성분 1의 전체 몰분율
A = np.array([13.8183, 13.8587])
B = np.array([2477.07, 2991.32])
C = np.array([233.21, 216.64])

# 2. 초기 추정값 설정 (x1, y1, P)
x0 = [0.1, 0.6, 50.0]

# 3. fsolve를 통한 비선형 방정식 풀이
# fsolve는 기본적으로 함수(args) 형식을 지원하므로 args로 나머지 파라미터를 넘깁니다.
solution = fsolve(binflash, x0, args=(v, T, z, A, B, C))

# 4. 결과 출력
x1_res, y1_res, P_res = solution

print(f"--- Flash Calculation Results (T={T}C) ---")
print(f"Liquid mole fraction (x1): {x1_res:.4f}")
print(f"Vapor mole fraction (y1): {y1_res:.4f}")
print(f"Pressure (P): {P_res:.2f}")