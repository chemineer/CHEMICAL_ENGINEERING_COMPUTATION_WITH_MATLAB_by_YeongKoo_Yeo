import numpy as np

def congrad(A, x0, b):
    """
    공액구배법을 사용하여 Ax = b를 풉니다.
    
    입력:
    A: 계수 행렬 (n x n)
    x0: 초기 추정값
    b: 상수 벡터
    """
    # 설정 및 초기화
    tol = 1e-6
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float).flatten()
    x = np.array(x0, dtype=float).flatten()
    n = len(b)
    
    # 1. 행렬식 확인 (Singularity 체크)
    if np.linalg.det(A) == 0:
        print(f"Rank = {np.linalg.matrix_rank(A)}")
        raise ValueError("Matrix A is singular.")
        
    # 2. 대각 우세(Diagonal Dominancy) 확인
    for k in range(n):
        diag_val = abs(A[k, k])
        off_diag_sum = np.sum(np.abs(A[k, :])) - diag_val
        if off_diag_sum > diag_val:
            print("Warning: A is not diagonally dominant.")
            # 원본 MATLAB 로직에 따라 대각 우세가 아니면 중단
            return None

    # 3. 반복 계산 (Iteration)
    # 초기 잔차(s)와 방향 벡터(d) 설정
    s = b - np.dot(A, x)
    d = np.copy(s)
    
    for k in range(n):
        v = np.dot(A, d)
        
        # 알파(step size) 계산: (d · s) / (d · v)
        alpha = np.dot(d, s) / np.dot(d, v)
        
        # 해 업데이트
        x = x + alpha * d
        
        # 새로운 잔차 계산
        s = b - np.dot(A, x)
        
        # 수렴 확인 (L2 Norm 사용)
        if np.sqrt(np.dot(s, s)) < tol:
            print(f"Converged at iteration: {k + 1}")
            return x
        else:
            # 베타 계산 및 방향 벡터 업데이트
            beta = -np.dot(s, v) / np.dot(d, v)
            d = s + beta * d
            
    return x

# --- 사용 예시 ---
if __name__ == "__main__":
    # 대칭이며 대각 우세인 행렬 예시
    A_mat = [[2, -1, 0, 0, 0, 0],
             [-1, 2, -1, 0, 0, 0],
             [0, -1, 2, -1, 0, 0],
             [0, 0, -1, 2, -1, 0],
             [0, 0, 0, -1, 2, -1],
             [0, 0, 0, 0, -1, 2]]
    b_vec = [0, 0, 0, 0, 0, 4]
    initial_guess = [0, 0, 0, 0, 0, 0]
    
    solution = congrad(A_mat, initial_guess, b_vec)
    if solution is not None:
        print(f"Solution: {solution}")