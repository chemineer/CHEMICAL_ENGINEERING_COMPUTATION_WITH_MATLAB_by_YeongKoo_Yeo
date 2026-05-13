import numpy as np
import scipy.linalg as la

def parabDbc(nx, nt, dx, dt, alpha, u0, bci, bcf):
    """
    1차원 포물형 편미분 방정식(PDE)을 상수 디리클레 경계 조건으로 해석합니다.
    """
    r = alpha * dt / dx**2
    nx = int(nx)
    nt = int(nt)
    
    # 행렬 A 초기화 (nx-1 x nx-1)
    A = np.zeros((nx - 1, nx - 1))
    u = np.zeros((nt + 1, nx + 1))
    
    # 경계 조건 및 초기 조건 설정
    u[:, 0] = bci
    u[:, nx] = bcf
    u[0, :] = u0
    
    # 행렬 A 구성 (크랭크-니콜슨 혹은 암시적 방법 계수)
    # 원본 코드 로직에 따라 삼중대각행렬 구성
    A[0, 0] = 1 + 2 * r
    if nx > 2:
        A[0, 1] = -r
        
    for i in range(1, nx - 2):
        A[i, i] = 1 + 2 * r
        A[i, i-1] = -r
        A[i, i+1] = -r
        
    if nx > 2:
        A[nx-2, nx-3] = -r
        A[nx-2, nx-2] = 1 + 2 * r
        
    # 벡터 b 초기화
    b = np.zeros(nx - 1)
    b[0] = u0[1] + u0[0] * r
    for i in range(1, nx - 2):
        b[i] = u0[i+1]
    b[nx-2] = u0[nx-1] + u0[nx] * r
    
    # LU 분해
    lu, piv = la.lu_factor(A)
    
    # 시간 단계 반복
    for j in range(1, nt + 1):
        x = la.lu_solve((lu, piv), b)
        u[j, 1:nx] = x
        
        # 다음 단계 b 업데이트
        b = x.copy()
        b[0] += bci * r
        b[nx-2] += bcf * r
        
    return u, r