import numpy as np

def gausselm(A, b):
    """
    가우스 소거법을 사용하여 Ax = b를 풉니다.
    
    입력:
    A: n x n 행렬 (리스트 또는 numpy array)
    b: n x 1 열 벡터 (리스트 또는 numpy array)
    """
    # 입력을 numpy 배열로 변환하고 부동소수점 타입(float) 확인
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float).flatten() # b를 1차원 배열로 처리
    n = len(b)
    
    # 1. Elimination phase (소거 단계)
    for k in range(n - 1):
        for i in range(k + 1, n):
            if A[i, k] != 0:
                # 피벗 계수 계산
                c = A[i, k] / A[k, k]
                
                # 행 연산: A[i, k+1:] = A[i, k+1:] - c * A[k, k+1:]
                A[i, k+1:] = A[i, k+1:] - c * A[k, k+1:]
                
                # 결과 벡터 b 업데이트
                b[i] = b[i] - c * b[k]
                
    # 2. Solution phase: back substitution (후진 대입법)
    x = np.zeros(n)
    for k in range(n - 1, -1, -1):
        # x[k] = (b[k] - (A[k, k+1:] * x[k+1:]의 합)) / A[k, k]
        x[k] = (b[k] - np.dot(A[k, k+1:], x[k+1:])) / A[k, k]
        
    return x

# --- 사용 예시 ---
if __name__ == "__main__":
    # Ax = b 시스템 설정
    A_mat = [[2, -1, 0, 0, 0, 0],
             [-1, 2, -1, 0, 0, 0],
             [0, -1, 2, -1, 0, 0],
             [0, 0, -1, 2, -1, 0],
             [0, 0, 0, -1, 2, -1],
             [0, 0, 0, 0, -1, 2]]
    b_vec = [0, 0, 0, 0, 0, 4]
    
    solution = gausselm(A_mat, b_vec)
    print(f"Solution: {solution}")