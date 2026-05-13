import numpy as np

def simplexlp(A, b, c, constr):
    """
    simplexlp.m: 2단계 심플렉스 알고리즘을 이용한 LP 최소화 문제 해결
    Minimize f(x) = c*x subject to Ax (constr) b, x >= 0
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float).reshape(-1, 1)
    c = np.array(c, dtype=float).flatten()
    
    m, n = A.shape
    n1 = n
    nleq = 0
    neq = 0
    ncomv = 0
    
    # c의 길이가 n보다 짧을 경우 0으로 채움
    if len(c) < n:
        c = np.pad(c, (0, n - len(c)), 'constant')
    
    # 제약 조건 타입에 따른 슬랙/잉여 변수 처리
    for j in range(m):
        temx = np.zeros((m, 1))
        temx[j] = 1
        if constr[j] == '<':  # <=
            A = np.hstack([A, temx])
            nleq += 1
        elif constr[j] == '>':  # >=
            A = np.hstack([A, -temx])
        else:  # =
            neq += 1
            
    lenA = A.shape[1]
    
    if nleq == m:
        # 모든 제약 조건이 <= 인 경우 (Phase 1 생략 가능 케이스)
        c_aug = np.pad(c, (0, lenA - len(c)), 'constant')
        A_final = np.vstack([A, c_aug])
        b_zero = np.vstack([b, [[0]]])
        A_final = np.hstack([A_final, b_zero])
        
        maux, A_final, z = compsim(A_final, list(range(n1, lenA)), 1, 1)
    else:
        # 2단계 심플렉스 (Phase 1: 인공 변수 도입)
        A_phase1 = np.hstack([A, np.eye(m), b])
        
        # 보조 목적 함수 (인공 변수의 합 최소화)
        if m > 1:
            w = -np.sum(A_phase1[:m, :lenA], axis=0)
        else:
            w = -A_phase1[0, :lenA]
            
        c_aug = np.pad(c, (0, A_phase1.shape[1] - 1 - len(c)), 'constant')
        
        # 테이블로 결합
        row_c = np.zeros((1, A_phase1.shape[1]))
        row_c[0, :len(c_aug)] = c_aug
        
        row_w = np.zeros((1, A_phase1.shape[1]))
        row_w[0, :lenA] = w
        row_w[0, -1] = -np.sum(b)
        
        A_tab = np.vstack([A_phase1[:m, :], row_c, row_w])
        
        maux = list(range(lenA, lenA + m))
        mv = list(maux)
        
        # Phase 1 실행
        maux, A_tab, z = compsim(A_tab, maux, 2, 1)
        
        nc = A_tab.shape[1]
        x_check = np.zeros(nc - 1)
        for i, idx in enumerate(maux):
            if idx < len(x_check):
                x_check[idx] = A_tab[i, -1]
        
        xm = x_check[mv]
        incomv = set(maux).intersection(set(mv))
        
        if np.any(np.abs(xm) > 1e-9):
            print("\nEmpty feasible region\n")
            return None, None
        else:
            if incomv:
                ncomv = 1
        
        # Phase 2를 위해 인공 변수 열 제거 및 목적 함수 재설정
        A_final = np.hstack([A_tab[:m+1, :lenA], A_tab[:m+1, -1:]])
        maux, A_final, z = compsim(A_final, maux, 1, 2)
    
    if np.isinf(z):
        return None, None
        
    m_f, n_f = A_final.shape
    xopt = np.zeros(n_f - 1)
    for i, idx in enumerate(maux):
        if idx < n_f - 1:
            xopt[idx] = A_final[i, -1]
            
    xopt = xopt[:n1]
    fopt = -A_final[-1, -1]
    
    # 다중 해 및 중복 제약 조건 체크
    t = np.where(np.abs(A_final[-1, :-1]) < 1e-9)[0]
    if len(t) > (m_f - 1):
        print('There are infinite solutions')
    if ncomv == 1:
        print('Redundant constraint(s).')
        
    return xopt, fopt

def compsim(A, maux, k, ph):
    """심플렉스 알고리즘의 메인 루프 (Bland's rule 적용)"""
    m, n = A.shape
    mi, col = Bland(A[m-1, :n-1])
    
    while col is not None and mi < -1e-12:
        t = A[:m-k, col]
        if np.all(t <= 0):
            print(f"\n Unbounded optimal solution\n")
            return maux, A, float('-inf')
            
        row, small = minrtest(A[:m-k, n-1], A[:m-k, col])
        
        if row is not None:
            # 피벗팅
            pivot_val = A[row, col]
            A[row, :] = A[row, :] / pivot_val
            maux[row] = col
            
            for i in range(m):
                if i != row:
                    A[i, :] = A[i, :] - A[i, col] * A[row, :]
            
            mi, col = Bland(A[m-1, :n-1])
        else:
            break
            
    return maux, A, A[m-1, n-1]

def Bland(D):
    """Bland's rule을 적용하여 피벗 열 선택"""
    ind = np.where(D < -1e-12)[0]
    if ind.size > 0:
        j = ind[0]
        return D[j], j
    return None, None

def minrtest(a, b):
    """Minimum ratio test를 수행하여 피벗 행 선택"""
    a = a.flatten()
    b = b.flatten()
    # 분모가 0보다 큰 경우에 대해서만 비율 계산
    eligible = np.where(b > 1e-12)[0]
    
    if eligible.size == 0:
        return None, None
        
    ratios = a[eligible] / b[eligible]
    idx = np.argmin(ratios)
    row = eligible[idx]
    return row, ratios[idx]

# --- 사용 예시 ---
if __name__ == "__main__":
    # 예시: Minimize f(x)=2x1 + 3x2 + 2x3 - x4 + x5
    # subject to 3x1-3x2+4x3+2x4-x5=0, x1+x2+x3+3x4+x5=2
    A = [[3, -3, 4, 2, -1], [1, 1, 1, 3, 1]]
    b = [0, 2]
    c = [2, 3, 2, -1, 1]
    constr = '=='
    
    xopt, fopt = simplexlp(A, b, c, constr)
    print(f"Optimal x: {xopt}")
    print(f"Optimal f: {fopt}")