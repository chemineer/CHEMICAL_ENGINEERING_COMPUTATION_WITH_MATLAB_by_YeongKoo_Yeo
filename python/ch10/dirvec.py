import numpy as np

def dirvec(delfun, n, nc, ne, crit, x, nca, A):
    """
    dirvec: 제약 조건이 있는 최적화 문제에서 탐색 방향을 결정하는 함수
    """
    df = np.array(delfun(x)).flatten()  # 그레이디언트 (1D 배열)
    
    while True:
        if nc == 0:
            d = -df
            dn = np.linalg.norm(d)
            if dn < crit:
                break
            d = d / dn
            Bm = np.array([0.0])
            return dn, d, Bm
        else:
            # Am 행렬 구성
            Am = np.zeros((nc, nc))
            for j in range(nc):
                ic = nca[j]
                for jk in range(nc):
                    jc = nca[jk]
                    # 행렬 곱: A(ic, :) dot A(jc, :)
                    Am[j, jk] = np.dot(A[ic], A[jc])
            
            # Bm 계산 (연립방정식 풀이)
            Bm_vec = np.zeros(nc)
            for j in range(nc):
                ic = nca[j]
                Bm_vec[j] = -np.dot(A[ic], df)
            
            # MATLAB의 inv(Am) * Bm 대신 np.linalg.solve 사용 (수치적으로 더 안정적)
            Bm = np.linalg.solve(Am, Bm_vec)
            
            # 탐색 방향 d 계산
            d = -df.copy()
            for k in range(nc):
                kn = nca[k]
                d -= A[kn] * Bm[k]
        
        dn = np.linalg.norm(d)
        
        # 중단 조건 및 제약 조건 처리
        if nc == ne and dn <= crit:
            break
        
        if dn <= crit:
            # Bm 중에서 최소값 찾기 (ne+1 인덱스부터)
            # 파이썬 인덱스는 0부터 시작하므로 ne 인덱스부터 탐색
            search_range = Bm[ne:]
            if len(search_range) > 0:
                Bmin = np.min(search_range)
                imin = ne + np.argmin(search_range)
                
                if Bmin >= 0:
                    break
                else:
                    # nca 리스트 내 요소 교환 (remove constraint)
                    nca[imin], nca[nc-1] = nca[nc-1], nca[imin]
                    nc -= 1
            else:
                break
        else:
            d = d / dn
            break
            
    return dn, d, Bm