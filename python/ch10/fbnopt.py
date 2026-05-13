import math

def fbnopt(objfun, a, b, n):
    """
    fbnopt: 1-dimensional Fibonacci search method
    :param objfun: 목적 함수 (람다식 또는 함수)
    :param a, b: 초기 구간 양 끝점
    :param n: 반복 횟수 (reduction count)
    :return: (xopt, fopt, fint) - 최적점, 최적점에서의 함수값, 최종 구간 크기
    """
    # 황금비와 관련된 상수 계산
    sqrt5 = math.sqrt(5)
    v = (sqrt5 - 1) / 2
    w = (1 - sqrt5) / (1 + sqrt5)
    
    x1 = float(a)
    x4 = float(b)
    
    # 초기 alpha 계산
    alpha = v * (1 - w**n) / (1 - w**(n + 1))
    x3 = alpha * x4 + (1 - alpha) * x1
    f3 = objfun(x3)
    
    for k in range(1, n):  # MATLAB의 k = 1:n-1 반복문 대응
        if k == n - 1:
            x2 = 0.01 * x1 + 0.99 * x3
        else:
            x2 = alpha * x1 + (1 - alpha) * x4
            
        f2 = objfun(x2)
        
        # 구간 축소 로직
        if f2 < f3:
            x4 = x3
            x3 = x2
            f3 = f2
        else:
            x1 = x4
            x4 = x2
            # f4 = f2 # MATLAB 원본에 명시된 할당이지만 이후 사용되지 않음
            
        # 다음 단계 alpha 업데이트
        alpha = v * (1 - w**(n - k)) / (1 - w**(n - k + 1))
        
    x = x3
    f = f3
    fint = abs(x1 - x4)
    
    return x, f, fint

# --- 사용 예시 ---
if __name__ == "__main__":
    # 목적 함수 정의: -870x + 102x^2 - 5x^3
    fobj = lambda x: -870*x + 102*x**2 - 5*x**3
    
    a = 0
    b = 50
    n = 20
    
    xopt, fopt, interval_size = fbnopt(fobj, a, b, n)
    
    print(f"Optimal point (x): {xopt:.6f}")
    print(f"Function value (f): {fopt:.6f}")
    print(f"Final interval size: {interval_size:.6f}")