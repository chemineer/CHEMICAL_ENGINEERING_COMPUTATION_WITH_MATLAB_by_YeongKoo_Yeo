import numpy as np

def barnslp(A, b, c, tol):
    """
    Barnes의 내점법을 사용한 선형 계획법 솔버
    """
    m, n_orig = A.shape
    
    # 1. c 벡터의 길이를 원래 변수 개수(n_orig)에 맞춤
    c_full = np.zeros(n_orig)
    c_full[:len(c)] = c
    
    # 2. 인공 변수를 위한 열(aplus1)과 큰 비용(cplus1) 추가
    aplus1 = (b.flatten() - np.sum(A, axis=1)).reshape(-1, 1)
    cplus1 = 1e6
    
    # 행렬 A와 벡터 c를 확장
    A_ext = np.hstack([A, aplus1])
    c_ext = np.append(c_full, cplus1)
    
    n = A_ext.shape[1] # 확장된 변수의 총 개수
    x = np.ones(n)
    alpha = 0.0001
    lambda_vec = np.zeros(m)
    
    # 루프 내부에서 사용할 변수명을 A_ext와 c_ext로 통일하거나 
    # 아래와 같이 A, c에 다시 할당합니다.
    A = A_ext
    c = c_ext
    
    # Main step
    while True:
        # 목적 함수 값과 제약 조건의 이중성 비교 (종료 조건)
        if abs(np.dot(c, x) - np.dot(lambda_vec, b)) <= tol:
            break
            
        x2 = x * x
        D = np.diag(x)
        D2 = np.diag(x2)
        AD2 = A @ D2
        
        # 람다 계산 (선형 방정식 풀이)
        lambda_vec = np.linalg.solve(AD2 @ A.T, AD2 @ c.T)
        
        # 이중 잔차 계산
        dualres = c - A.T @ lambda_vec
        normres = np.linalg.norm(D @ dualres)
        
        # 비율 계산
        ratio = np.full(n, np.inf)
        for i in range(n):
            denom = x[i] * (c[i] - np.dot(A[:, i], lambda_vec))
            if dualres[i] > 0:
                ratio[i] = normres / denom
                
        # 업데이트 단계
        R = np.min(ratio) - alpha
        x = x - R * D2 @ dualres / normres
        
        # 기저 변수(Basic variables) 확인
        basic = np.where(x > tol)[0]
        basiscount = len(basic)
        
        # 비퇴화 문제일 경우 솔루션 도출
        if basiscount == m:
            B = A[:, basic]
            # primalsol = b'/B'와 같은 의미 (B @ x = b 해결)
            xopt_full = np.linalg.lstsq(B, b, rcond=None)[0]
            
            # 최종 결과 정리
            xopt = xopt_full
            fopt = np.dot(c[basic], xopt)
            return xopt, fopt, basic + 1 # 1-based index 반환
            
    # 루프 탈출 후 최종 결과
    xopt = x[basic]
    fopt = np.dot(c, x)
    return xopt, fopt, basic + 1