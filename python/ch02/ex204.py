import numpy as np
from scipy.sparse.linalg import cg

# 1. 행렬 A와 벡터 b 정의
A = np.array([
    [ 4, -1, -1,  0],
    [-1,  4,  0, -1],
    [-1,  0,  4, -1],
    [ 0, -1, -1,  4]
], dtype=float)

b = np.array([45, 35, 55, 45], dtype=float)

# 2. 공액구배법(CG) 수행
# scipy.sparse.linalg.cg(A, b)는 MATLAB의 pcg(A, b)와 대응됩니다.
# x: 해 벡터, info: 수렴 정보 (0이면 성공)
x, info = cg(A, b)

# 3. 결과 출력
if info == 0:
    print("Optimization terminated successfully.")
    print("Solution x:")
    print(x)
else:
    print(f"Convergence failed. Info code: {info}")

# --- 참고: 직접 역행렬과 비교 (검증용) ---
# print("\nVerification (A \ b):")
# print(np.linalg.solve(A, b))