import numpy as np

def newtintp(x, y, xi):
    """
    뉴턴 분점차 보간법(Newton Interpolation) 구현 함수
    
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
    
    # 1. 뉴턴 다항식의 계수(분점차) 결정
    # df 테이블 생성: df[행, 열]
    df = np.zeros((n, n))
    a = np.zeros(n)
    
    a[0] = y[0]
    
    # 첫 번째 차수 분점차 계산 (MATLAB의 df(k, 1))
    for k in range(n - 1):
        df[k, 0] = (y[k+1] - y[k]) / (x[k+1] - x[k])
        
    # 고차 분점차 계산 (MATLAB의 df(k, j))
    for j in range(1, n - 1): # 열 인덱스 (1부터 n-2까지)
        for k in range(n - (j + 1)): # 행 인덱스
            df[k, j] = (df[k+1, j-1] - df[k, j-1]) / (x[k+j+1] - x[k])
            
    # 계수 a 추출 (첫 번째 행의 값들)
    for k in range(1, n):
        a[k] = df[0, k-1]
        
    # 2. 뉴턴 다항식을 이용한 보간 수행
    yi = np.zeros(m)
    for k in range(m):
        s = np.zeros(n)
        s[0] = 1.0
        yi[k] = a[0]
        for j in range(1, n):
            s[j] = (xi[k] - x[j-1]) * s[j-1]
            yi[k] += a[j] * s[j]
            
    return yi

# --- 테스트 코드 ---
if __name__ == "__main__":
    # 데이터 예시
    x_data = [1, 2, 4]
    y_data = [1, 4, 16]
    x_target = [1.5, 3.0]
    
    result = newtintp(x_data, y_data, x_target)
    print(f"보간 위치: {x_target}")
    print(f"보간 결과: {result}")