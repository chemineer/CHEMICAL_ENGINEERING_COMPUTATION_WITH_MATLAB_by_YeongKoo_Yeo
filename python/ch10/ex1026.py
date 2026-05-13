import numpy as np
from scipy.optimize import minimize

def main():
    # 1. 목적 함수 정의 (MATLAB: f = @(x) 2*x^2*sin(x) + exp(-x))
    # 단일 변수 x에 대한 함수입니다.
    def f(x):
        return 2 * (x**2) * np.sin(x) + np.exp(-x)

    # 2. 초기값 설정 (MATLAB: x0 = 1)
    x0 = 1.0

    # 3. fminsearch에 대응하는 최적화 (Nelder-Mead 알고리즘)
    # 미분값을 사용하지 않는 심플렉스 방식입니다.
    res_search = minimize(f, x0, method='Nelder-Mead')
    
    # 4. fminunc에 대응하는 최적화 (BFGS 등 경사하강 기반 알고리즘)
    # 일반적으로 SciPy의 기본값이나 BFGS가 fminunc와 유사하게 작동합니다.
    res_unc = minimize(f, x0, method='BFGS')

    # 5. 결과 출력
    print("=== SciPy Optimization Results (Ex 10.26) ===")
    
    print("\n[Method: Nelder-Mead (fminsearch 대체)]")
    print(f"최적해 x: {res_search.x[0]:.6f}")
    print(f"최적 함수값 f(x): {res_search.fun:.6f}")
    
    print("\n[Method: BFGS (fminunc 대체)]")
    print(f"최적해 x: {res_unc.x[0]:.6f}")
    print(f"최적 함수값 f(x): {res_unc.fun:.6f}")

if __name__ == "__main__":
    main()