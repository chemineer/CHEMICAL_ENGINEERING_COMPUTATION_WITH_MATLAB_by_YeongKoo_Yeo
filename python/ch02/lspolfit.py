import numpy as np

def lspolfit(x, y, m):
    """
    최소제곱 다항식 적합(Least-squares polynomial fitting)을 수행합니다.
    
    입력:
    x: 독립 변수 벡터
    y: 측정된 종속 변수 벡터
    m: 다항식의 차수
    
    출력:
    p: m차 다항식의 계수 벡터 (내림차순)
    """
    x = np.array(x)
    y = np.array(y)
    n = len(x)
    
    if n != len(y):
        raise ValueError('x와 y의 길이는 같아야 합니다.')
        
    # 다항식 차수 m에 대해 행렬 크기는 (m+1, m+1)
    M = m + 1
    A = np.zeros((M, M))
    
    # x의 거듭제곱 합 계산 (s_k = sum(x^k))
    # 최대 2*m 차수까지 필요함
    s = np.zeros(2 * m)
    for k in range(2 * m):
        s[k] = np.sum(x**(k + 1))
        
    # 행렬 A 정의 (Normal Equation Matrix)
    A[0, 0] = n
    A[0, 1:M] = s[0:M-1]
    for k in range(1, M):
        # MATLAB의 k번째 행은 파이썬의 k-1 인덱스
        A[k, :] = s[k-1 : k-1+M]
        
    # 벡터 b 정의 (sum(x^k * y))
    b = np.zeros(M)
    for k in range(M):
        b[k] = np.sum((x**k) * y)
        
    # 선형 시스템 Ac = b 풀기 (MATLAB의 c = A\b)
    c = np.linalg.solve(A, b)
    
    # 계수를 내림차순으로 정렬 (MATLAB의 p = c(end:-1:1))
    p = c[::-1]
    
    return p
    
