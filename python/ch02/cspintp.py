import numpy as np

def cspintp(x, y, xi):
    """
    Cubic Spline Interpolation (삼차 스플라인 보간법) 구현
    
    입력:
    x: 독립 변수 데이터 (알려진 점들)
    y: 종속 변수 데이터
    xi: 보간 값을 구하고자 하는 목표 지점들
    
    출력:
    yi: xi 지점들에서의 보간 결과값
    """
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)
    xi = np.array(xi, dtype=float)
    
    n = len(x)
    m = len(xi)
    
    if n != len(y):
        raise ValueError("x와 y의 길이는 같아야 합니다.")
    
    # 1. Tridiagonal System 구성을 위한 파라미터 계산
    h = np.zeros(n - 1)
    w = np.zeros(n - 1)
    for k in range(n - 1):
        h[k] = x[k+1] - x[k]
        w[k] = (y[k+1] - y[k]) / h[k]
        
    v = np.zeros(n - 2)
    d = np.zeros(n - 2)
    for k in range(n - 2):
        v[k] = w[k+1] - w[k]
        d[k] = 2 * (h[k] + h[k+1])
        
    U = np.zeros(n - 2)
    L = np.zeros(n - 2)
    for k in range(n - 3):
        U[k] = h[k+1]
        L[k+1] = U[k]
    
    L[0] = 0
    U[n-3] = 0
    
    # 2. Tridiagonal System 해결 (계수 a 계산)
    v_proc = v.copy()
    U_proc = U.copy()
    
    v_proc[0] = v_proc[0] / d[0]
    U_proc[0] = U_proc[0] / d[0]
    
    for k in range(1, n - 3):
        dn = d[k] - L[k] * U_proc[k-1]
        U_proc[k] = U_proc[k] / dn
        v_proc[k] = (v_proc[k] - L[k] * v_proc[k-1]) / dn
        
    # 마지막 행 처리
    last_idx = n - 3
    v_proc[last_idx] = (v_proc[last_idx] - L[last_idx] * v_proc[last_idx-1]) / \
                       (d[last_idx] - L[last_idx] * U_proc[last_idx-1])
    
    a_coeffs = np.zeros(n - 2)
    a_coeffs[last_idx] = v_proc[last_idx]
    for k in range(n - 4, -1, -1):
        a_coeffs[k] = v_proc[k] - U_proc[k] * a_coeffs[k+1]
        
    # 3. 보간 식을 위한 계수 b, c 계산
    b = np.zeros(n - 1)
    c = np.zeros(n - 1)
    
    b[0] = y[0] / h[0]
    c[0] = y[1] / h[0] - a_coeffs[0] * h[0]
    
    for k in range(1, n - 2):
        b[k] = y[k] / h[k] - a_coeffs[k-1] * h[k]
        c[k] = y[k+1] / h[k] - a_coeffs[k] * h[k]
        
    b[n-2] = y[n-2] / h[n-2] - a_coeffs[n-3] * h[n-2]
    c[n-2] = y[n-1] / h[n-2]
    
    # a 벡터의 맨 앞에 0 추가 (Natural Spline 조건 반영)
    full_a = np.insert(a_coeffs, 0, 0)
    
    # 4. xi 지점들에서 보간 수행
    s = np.zeros(m)
    for k in range(m):
        # 현재 xi가 어느 구간 [x(id), x(id+1)]에 속하는지 확인
        if xi[k] > x[-1] or xi[k] < x[0]:
            raise ValueError(f"xi={xi[k]}는 보간 범위를 벗어났습니다.")
            
        id_idx = 0
        for j in range(n - 1):
            if xi[k] >= x[j]:
                id_idx = j
        
        h_val = x[id_idx+1] - x[id_idx]
        
        if id_idx == 0:
            s[k] = full_a[1] * (xi[k] - x[0])**3 / h_val + \
                   b[0] * (x[1] - xi[k]) + c[0] * (xi[k] - x[0])
        elif id_idx == n - 2:
            s[k] = full_a[n-2] * (x[n-1] - xi[k])**3 / h_val + \
                   b[n-2] * (x[n-1] - xi[k]) + c[n-2] * (xi[k] - x[n-2])
        else:
            term1 = (full_a[id_idx] * (x[id_idx+1] - xi[k])**3 + \
                     full_a[id_idx+1] * (xi[k] - x[id_idx])**3) / h_val
            term2 = b[id_idx] * (x[id_idx+1] - xi[k]) + \
                    c[id_idx] * (xi[k] - x[id_idx])
            s[k] = term1 + term2
            
    return s

# --- 테스트 실행 ---
if __name__ == "__main__":
    x_pts = [0, 1, 2, 3]
    y_pts = [0, 1, 0, 1]
    xi_pts = [0.5, 1.5, 2.5]
    
    yi_pts = cspintp(x_pts, y_pts, xi_pts)
    print("보간 결과:", yi_pts)