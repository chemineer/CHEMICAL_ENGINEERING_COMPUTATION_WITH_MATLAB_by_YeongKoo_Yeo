import numpy as np

def GaussSeidel(A, b, rho):
    """
    가우스-사이델 반복법을 사용하여 Ax = b를 풉니다.
    
    입력:
    A: 계수 행렬 (n x n)
    b: 상수 벡터
    rho: 이완 계수 (Relaxation factor)
    """
    # 초기 설정
    itmax = 500
    tol = 1e-8
    
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float).flatten()
    n = len(b)
    nr, nc = A.shape
    
    # 예외 처리: 정방 행렬 및 크기 일치 확인
    if nr != nc:
        raise ValueError("Matrix A is not square.")
    if nr != n:
        raise ValueError("Matrix A and vector b are not consistent.")
    
    # 행렬식 확인 (Singularity 체크)
    if np.linalg.det(A) == 0:
        print(f"Rank = {np.linalg.matrix_rank(A)}")
        raise ValueError("A is singular.")
        
    # 대각 우세(Diagonally Dominant) 확인 및 행렬 M 구성
    M = np.copy(A)
    x = np.zeros(n)
    s = np.zeros(n)
    
    for k in range(n):
        diag_val = abs(A[k, k])
        off_diag_sum = np.sum(np.abs(A[k, :])) - diag_val
        if off_diag_sum > diag_val:
            print("Warning: A is not diagonally dominant.")
            # MATLAB 코드에서는 여기서 return하지만, 수렴 가능성이 있으므로 계속 진행하도록 둠
            
        # 행렬 M과 벡터 s 정규화 (A(k,k)로 나누기)
        M[k, :] = A[k, :] / A[k, k]
        M[k, k] = 0 # 자기 자신(대각 성분)은 계산에서 제외하도록 0으로 설정
        s[k] = b[k] / A[k, k]

    iter_count = 0
    err = np.zeros(n)
    
    # 반복 계산 시작
    while True:
        x0 = np.copy(x)
        
        for k in range(n):
            # 가우스-사이델 핵심: 현재까지 업데이트된 x값을 바로 사용
            # x(k) = rho * (s(k) - sum(M(k,j) * x(j))) + (1-rho) * x0(k)
            x[k] = rho * (s[k] - np.dot(M[k, :], x)) + (1 - rho) * x0[k]
            
            if x[k] != 0:
                err[k] = abs((x[k] - x0[k]) / x[k])
        
        iter_count += 1
        
        # 종료 조건: 허용 오차 도달 또는 최대 반복 횟수 초과
        if np.max(err) <= tol or iter_count >= itmax:
            break
            
    print(f"Number of iterations: {iter_count}")
    return x

# --- 사용 예시 ---
if __name__ == "__main__":
    A_mat = [[2, -1, 0, 0, 0, 0],
             [-1, 2, -1, 0, 0, 0],
             [0, -1, 2, -1, 0, 0],
             [0, 0, -1, 2, -1, 0],
             [0, 0, 0, -1, 2, -1],
             [0, 0, 0, 0, -1, 2]]
    b_vec = [0, 0, 0, 0, 0, 4]
    rho_val = 0.8 # 1이면 기본 Gauss-Seidel
    
    result = GaussSeidel(A_mat, b_vec, rho_val)
    print(f"Solution: {result}")