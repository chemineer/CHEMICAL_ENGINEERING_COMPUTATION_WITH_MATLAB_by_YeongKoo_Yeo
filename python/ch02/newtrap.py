def newtrap(f, df, x0):
    """
    뉴턴-랩슨법을 사용하여 방정식의 근을 찾습니다.
    
    입력:
    f: f(x)를 반환하는 함수 핸들
    df: f'(x) (도함수)를 반환하는 함수 핸들
    x0: 초기 추정값
    
    출력:
    x: f(x)의 제로점(근)
    """
    tol = 1e-8
    kmax = 10000
    
    # 초기값 계산
    f0 = f(x0)
    df0 = df(x0)
    
    # 초기값이 근인 경우 처리
    if f0 == 0:
        return x0
    
    k = 0
    x1 = x0 # x1 초기화
    
    for i in range(1, kmax + 1):
        k = i
        
        # 뉴턴-랩슨 공식: x1 = x0 - f(x0)/f'(x0)
        if df0 == 0:
            print("Derivative is zero. No solution found.")
            return None
            
        x1 = x0 - f0 / df0
        f1 = f(x1)
        df1 = df(x1)
        
        # 수렴 조건 확인 (함수값이 tol 이하이거나, x의 변화량이 tol 이하일 때)
        if abs(f1) <= tol or abs(x1 - x0) <= tol:
            break
            
        if k >= kmax:
            print("The Newton-Raphson method not converged.")
            break
            
        # 다음 반복을 위해 값 업데이트
        x0 = x1
        f0 = f1
        df0 = df1
            
    print(f"Number of iterations by Newton-Raphson: {k}")
    return x1

# --- 사용 예시 ---
if __name__ == "__main__":

    my_f = lambda r: 9.496e3 * (1 - 12*r**2 + 16*r**3) - 1800
    my_df = lambda r: 9.496e3 * (-24*r + 48*r**2)
    
    root = newtrap(my_f, my_df, 0.3)
    print(f"Result: {root}")