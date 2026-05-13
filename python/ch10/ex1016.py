import numpy as np
from fun import fun
from dfun import dfun
from sqpopt import sqpopt

def main():
    # 1. 초기값 설정 (MATLAB: x0 = [2 3]', lam0 = 1, mu0 = 3)
    # x0: 초기 결정 변수 [x1, x2]
    x0 = np.array([2, 3], dtype=float)
    
    # lam0: 등식 제약 조건의 개수
    lam0 = 1
    
    # mu0: 부등식 제약 조건의 개수
    # 참고: fun.py에는 g1~g4까지 4개의 부등식 제약 조건이 정의되어 있습니다.
    # MATLAB 코드 예시에 따라 3으로 설정하거나, 실제 정의된 개수인 4로 설정할 수 있습니다.
    mu0 = 4 
    
    # crit: 수렴 판정 기준
    crit = 1e-6
    
    try:
        # 2. SQP 최적화 실행
        # xopt: 최적해, fopt: 최적 목적 함수 값, iter_count: 반복 횟수
        xopt, fopt, iter_count = sqpopt(fun, dfun, x0, lam0, mu0, crit)
        
        # 3. 결과 출력
        print("=== SQP Optimization Result ===")
        print(f"최적해 (xopt): {xopt}")
        print(f"최적 목적 함수 값 (fopt): {fopt:.6f}")
        print(f"반복 횟수 (iterations): {iter_count}")
        
        # 4. 제약 조건 검증
        fv = fun(xopt)
        print("\n--- Verification ---")
        print(f"목적 함수 값: {fv[0]:.6f}")
        print(f"등식 제약 조건 h(x)=0: {fv[1]:.6e}")
        print(f"부등식 제약 조건 g(x)>=0: {fv[2:]}")

    except Exception as e:
        print(f"최적화 중 오류가 발생했습니다: {e}")

if __name__ == "__main__":
    main()