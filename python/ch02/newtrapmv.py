import numpy as np

def newtrapmv(fun, x0, tol=1e-8, kmax=1000):
    """
    다변수 방정식 시스템의 해를 구하기 위한 뉴턴-랩슨법 구현
    
    입력:
    - fun: f(x)를 반환하는 함수 (리스트나 numpy 배열 반환)
    - x0: 초기 추측값 (리스트나 numpy 배열)
    
    출력:
    - x: 방정식의 해
    """
    x = np.array(x0, dtype=float).flatten()
    
    for k in range(1, kmax + 1):
        # 자코비언 행렬 J와 함수값 f 계산
        J, f = jacob(fun, x)
        
        # 1. 수렴 조건 확인 (RMS 오차)
        if np.sqrt(np.dot(f, f) / len(x)) < tol:
            print(f"Number of iterations: {k}")
            return x
        
        # 2. 증분 dx 계산 (J * dx = -f 해결)
        try:
            dx = np.linalg.solve(J, -f)
        except np.linalg.LinAlgError:
            print("Jacobian is singular. Newton-Raphson cannot proceed.")
            return None
            
        # x 업데이트
        x = x + dx
        
        # 3. 상대적 수렴 조건 확인
        if np.sqrt(np.dot(dx, dx) / len(x)) < tol * max(np.max(np.abs(x)), 1.0):
            print(f"Number of iterations: {k}")
            return x
            
    print('The Newton-Raphson method does not converge.')
    return x

def jacob(fun, x):
    """
    자코비언 행렬 J와 f(x)를 수치적으로 계산
    """
    hx = 1e-4
    n = len(x)
    f0 = np.array(fun(x), dtype=float)
    J = np.zeros((n, n))
    
    for k in range(n):
        tempx = x[k]
        x_plus_h = np.copy(x)
        x_plus_h[k] = tempx + hx
        
        f1 = np.array(fun(x_plus_h), dtype=float)
        J[:, k] = (f1 - f0) / hx
        
    return J, f0

# --- 사용 예시 ---
if __name__ == "__main__":
    # 예시 방정식: f1 = x^2 + y^2 - 4, f2 = e^x + y - 1
    def example_fun(x):
        f1 = np.cos(x[0]) + x[1]**2 + np.log(x[2]) - 8
        f2 = 4*x[0] + 3**x[1] - x[2]**3 + 2
        f3 = x[0] + x[1] + x[2] - 6
        return np.array([f1, f2, f3])
        

    initial_guess = np.array([1, 1, 1])
    solution = newtrapmv(example_fun, initial_guess)
    print(f"Solution: {solution}")