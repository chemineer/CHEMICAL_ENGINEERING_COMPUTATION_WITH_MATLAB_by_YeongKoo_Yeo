import numpy as np

def simps(f, a, b, n):
    """
    심슨 1/3 공식을 사용하여 f(x)를 a부터 b까지 적분합니다.
    
    입력:
    f: 적분할 함수 (함수 핸들)
    a, b: 적분 구간의 시작과 끝
    n: 격자점의 개수 (x1, ..., xn)
    
    출력:
    z: 적분 결과값
    """
    # 1. n이 짝수이면 홀수로 만듦 (심슨 공식은 구간 수가 짝수, 즉 점의 수가 홀수여야 함)
    if n % 2 == 0:
        n = n + 1
        
    h = (b - a) / (n - 1)
    s = f(a)
    
    # 2. 예외 처리 및 간단한 경우 계산
    if n <= 2:
        print("Too few subintervals.")
        return None
    
    if n == 3:
        z = h * (s + 4 * f((a + b) / 2) + f(b)) / 3
        return z
    
    # 3. 중간항 계산 (4 * 홀수 번째 인덱스 항들)
    # MATLAB: k = 2:2:n-1 -> Python: 인덱스 1, 3, 5... (0부터 시작하므로)
    for k in range(2, n, 2):
        x = a + h * (k - 1)
        s = s + 4 * f(x)
        
    # 4. 중간항 계산 (2 * 짝수 번째 인덱스 항들)
    # MATLAB: k = 3:2:n-2 -> Python: 인덱스 2, 4, 6...
    for k in range(3, n - 1, 2):
        x = a + h * (k - 1)
        s = s + 2 * f(x)
        
    s = s + f(b)
    z = h * s / 3
    
    return z