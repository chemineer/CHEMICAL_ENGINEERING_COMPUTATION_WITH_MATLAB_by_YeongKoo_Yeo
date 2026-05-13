def secant(fun, x1, x2):
    """
    할선법을 사용하여 방정식의 근을 찾습니다.
    
    입력:
    fun: f(x)를 반환하는 함수
    x1, x2: 근 근처의 초기 두 지점
    """
    tol = 1e-8 # 허용 오차
    kmax = 10000 # 최대 반복 횟수
    
    f1 = fun(x1) #
    f2 = fun(x2) #
    
    if f1 == 0: return x1 #
    if f2 == 0: return x2 #
    
    x3 = x2 # x3 초기화
    k = 0
    for i in range(1, kmax + 1):
        k = i
        # 할선법 공식: 다음 근사값 x3 계산
        x3 = x2 - f2 * (x2 - x1) / (f2 - f1)
        f3 = fun(x3) #
        
        if abs(f3) <= tol: # 수렴 조건 확인
            break
            
        if k >= kmax: # 최대 반복 도달 시 알림
            print('The secant method not converged.')
            break
            
        # 다음 반복을 위한 값 업데이트
        x1 = x2
        x2 = x3
        f1 = f2
        f2 = f3
        
    print(f'Number of iterations by secant: {k}') #
    return x3
    
# --- 사용 예시 ---
if __name__ == "__main__":

    my_fun = lambda r: 9.496e3 * (1 - 12*r**2 + 16*r**3) - 1800
    root = secant(my_fun, 0, 0.5)
    print(f"Result: {root}")