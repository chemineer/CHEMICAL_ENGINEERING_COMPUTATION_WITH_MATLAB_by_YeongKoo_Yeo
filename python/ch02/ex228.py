import numpy as np

# 1. 데이터 정의 (온도 T 및 결과값 b)
T = np.array([313.15, 363.15, 413.15])
b = np.array([8.2583, 7.4799, 6.9284])

# 2. 행렬 A 구성 (T^2, T, 1)
# MATLAB: A = [T.^2 T ones(3,1)]
A = np.column_stack([T**2, T, np.ones(3)])

# 3. 선형 방정식 Ax = b 풀기 (계수 x 계산)
# MATLAB: x = A\b
x = np.linalg.solve(A, b)

# 4. 결과 출력 (다항식 계수)
print("다항식 계수 (x1, x2, x3):")
print(x)

# 5. 373.15K에서의 보간값 계산
# MATLAB: s = x(1)*373.15^2 + x(2)*373.15 + x(3)
T_target = 373.15
s = x[0] * T_target**2 + x[1] * T_target + x[2]

print(f"\n온도 {T_target}K에서의 보간된 값 (s):")
print(f"{s:.10f}")