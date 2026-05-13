import numpy as np

# 1. 행렬 A 정의 (4x3 행렬)
A = np.array([
    [-2, -1,  3],
    [ 4, -5,  7],
    [ 6,  2, -5],
    [-3,  2,  1]
], dtype=float)

# 2. 의사역행렬(Pseudo-Inverse) 계산
Ap = np.linalg.pinv(A)

# 3. 결과 출력
print("Matrix A:")
print(A)

print("\nPseudo-Inverse of A (Ap):")
print(Ap)

# --- 검증 (A * Ap * A = A 성질 확인) ---
# print("\nVerification (A * Ap * A):")
# print(np.allclose(A, A @ Ap @ A))