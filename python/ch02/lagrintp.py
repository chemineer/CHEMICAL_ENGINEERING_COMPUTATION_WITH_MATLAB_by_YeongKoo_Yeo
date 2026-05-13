import numpy as np

def lagrintp(x, y, xi):
    """
    라그랑주 보간법을 구현한 함수
    
    입력:
    x: 독립 변수 벡터 (데이터 포인트)
    y: 종속 변수 벡터 (데이터 포인트)
    xi: 보간을 수행할 독립 변수 지점들
    
    출력:
    yi: xi 지점들에서의 보간된 값들
    """
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)
    xi = np.array(xi, dtype=float)
    
    n = len(x)
    m = len(xi)
    
    if n != len(y):
        raise ValueError("x와 y의 길이는 같아야 합니다.")
        
    # 라그랑주 다항식의 계수 a 계산
    a = np.zeros(n)
    for k in range(n):
        dx_k = 1.0
        for j in range(n):
            if j != k:
                dx_k *= (x[k] - x[j])
        a[k] = y[k] / dx_k
        
    # xi 지점들에서 보간 수행
    yi = np.zeros(m)
    for i in range(m):
        for j in range(n):
            q_j = 1.0
            for k in range(n):
                if j != k:
                    q_j *= (xi[i] - x[k])
            yi[i] += a[j] * q_j
            
    return yi

# --- 사용 예시 ---
if __name__ == "__main__":
    # 데이터 포인트 (예: y = x^2)
    x_data = [1, 2, 4]
    y_data = [1, 4, 16]
    
    # 보간할 지점
    x_interp = [1.5, 3.0]
    
    # 함수 실행
    y_interp = lagrintp(x_data, y_data, x_interp)
    
    print(f"보간 지점: {x_interp}")
    print(f"결과 값: {y_interp}")