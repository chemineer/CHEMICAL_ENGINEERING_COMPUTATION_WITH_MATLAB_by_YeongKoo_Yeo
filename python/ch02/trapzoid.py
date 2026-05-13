import numpy as np

def trapzoid(f, a, b, n):
    """
    함수 f를 a부터 b까지 n개의 구간으로 나누어 적분합니다.
    """
    h = (b - a) / n
    s = f(a)
    
    for k in range(1, n):
        x = a + h * k
        s += 2 * f(x)
        
    s += f(b)
    z = h * s / 2
    return z

def trapzoidat(x, y):
    """
    주어진 데이터 포인트 (x, y) 집합을 사다리꼴 공식으로 적분합니다.
    """
    n = len(x)
    s = 0.0
    
    for k in range(n - 1):
        # 각 구간의 너비와 높이의 평균을 곱하여 더함
        s += (y[k] + y[k+1]) * (x[k+1] - x[k]) / 2
        
    return s

# --- 사용 예시 ---

# 1. 함수 적분 예시: f(x) = x^2 을 0부터 1까지 적분
func = lambda x: x**2
result_func = trapzoid(func, 0, 1, 100)
print(f"함수 적분 결과: {result_func}")

# 2. 데이터 적분 예시
x_data = np.linspace(0, 1, 11) # 0부터 1까지 11개 포인트
y_data = x_data**2            # y = x^2
result_data = trapzoidat(x_data, y_data)
print(f"데이터 적분 결과: {result_data}")