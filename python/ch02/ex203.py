import numpy as np

A = np.array([[4, 2, -1], 
              [-3, 1, 2], 
              [2, -4, 1]], dtype=float)
b = np.array([8, -6, 12], dtype=float)

# MATLAB의 x = A\b와 동일
x = np.linalg.solve(A, b)

print("Solution x:")
print(x)

# rcond=None은 최신 버전에서 권장되는 설정입니다.
x, residuals, rank, s = np.linalg.lstsq(A, b, rcond=None)
print(x)

# MATLAB의 x = inv(A)*b와 동일
x = np.linalg.inv(A) @ b
print(x)