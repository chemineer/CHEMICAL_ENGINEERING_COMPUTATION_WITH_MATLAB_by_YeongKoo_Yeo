import math

def bisectn(fun, x1, x2):
    """
    이분법을 사용하여 함수의 근을 찾습니다.
    
    입력:
    fun: f(x)를 반환하는 함수 (lambda 또는 정의된 함수)
    x1, x2: 근을 포함하는 구간의 양 끝값
    
    출력:
    x: f(x)의 제로점(근)
    """
    tol = 1e-8
    
    # 반복 횟수 계산: kmax = ceil(log2(구간길이 / 허용오차))
    kmax = math.ceil(math.log(abs(x2 - x1) / tol) / math.log(2))
    
    f1 = fun(x1)
    f2 = fun(x2)
    
    # 초기 끝값이 근인 경우 처리
    if f1 == 0:
        return x1
    if f2 == 0:
        return x2
    
    # 사이값 정리 확인 (f1과 f2의 부호가 같으면 에러)
    if f1 * f2 > 0:
        raise ValueError("The root is not located in [x1, x2].")
    
    k = 0 # 반복 횟수 추적용
    for i in range(1, kmax + 1):
        k = i
        x3 = (x1 + x2) / 2
        f3 = fun(x3)
        
        # 허용 오차 내에 도달하면 중단
        if abs(f3) <= tol:
            x = x3
            break
            
        # 부호 변화가 있는 구간으로 좁히기
        if f2 * f3 < 0:
            x1 = x3
            f1 = f3
        else:
            x2 = x3
            f2 = f3
            
    x = (x1 + x2) / 2
    print(f"Number of iterations by bisect: {k}")
    return x

# --- 사용 예시 ---
if __name__ == "__main__":

    my_fun = lambda r: 9.496e3 * (1 - 12*r**2 + 16*r**3) - 1800
    root = bisectn(my_fun, 0, 0.5)
    print(f"Result: {root}")