import numpy as np
from grgfun import grgfun
from grgopt import grgopt

def delgrgf(x):
    """목적 함수의 기울기(Gradient) 계산"""
    x1, x2, x3, x4 = x[0], x[1], x[2], x[3]
    df = np.zeros(len(x))
    df[0] = 2.3 * x1 - 4
    df[1] = 2 * x2 - 6
    df[2] = 4.6 * x3 - 20
    df[3] = -2.4 * x4 + 6
    return df

def delgrgg(x):
    """제약 조건의 야코비안(Jacobian) 행렬 계산"""
    x1, x2, x3, x4 = x[0], x[1], x[2], x[3]
    # 3개의 제약 조건과 7개의 변수에 대한 미분 행렬
    dg = np.array([
        [2*x1 + 1, 2*x2 - 1, 3*x3 + 1, 2*x4 - 1, 1, 0, 0],
        [2*x1 - 1, 4*x2,     2*x3,     4*x4 - 1, 0, 1, 0],
        [4*x1 + 2, 2*x2 - 1, 2*x3,     -1,       0, 0, 1]
    ])
    return dg

def main():
    # 초기값 설정 (x1~x4는 0, x5~x7은 슬랙 변수 초기값)
    x0 = np.array([0, 0, 0, 0, 7, 11, 6], dtype=float)
    
    # 변수의 하한(Lower bound) 및 상한(Upper bound)
    xl = np.array([-100, -100, -100, -100, 0, 0, 0], dtype=float)
    xu = 100 * np.ones(7)
    
    # 제어 파라미터
    kmax = 1000      # 최대 반복 횟수
    crit = 1e-4      # 수렴 판정 기준
    
    # GRG 최적화 실행
    xopt, fopt, iter_count = grgopt(grgfun, delgrgf, delgrgg, x0, xl, xu, kmax, crit)
    
    # 결과 출력
    print("=== GRG Optimization Result ===")
    print(f"최적해 (xopt): {xopt}")
    print(f"최적 목적 함수 값 (fopt): {fopt:.6f}")
    print(f"반복 횟수 (iterations): {iter_count}")
    
    # 제약 조건 위반 여부 확인 (검증용)
    _, g_vals = grgfun(xopt)
    print(f"최적점에서의 제약 조건 (g): {g_vals}")

if __name__ == "__main__":
    main()