import numpy as np
import math
from saopt import saopt

def main():
    # 1. 목적 함수 정의 (MATLAB: @fcy)
    # f = -cos(5*sqrt(sum(x.^2))) + 0.1*sum(x.^2) 수식을 구현합니다.
    def fcy(x):
        # x는 [x0, x1] 형태의 배열
        r2 = np.sum(x**2)
        term1 = -math.cos(5 * math.sqrt(r2))
        term2 = 0.1 * r2
        return term1 + term2

    # 2. 초기 매개변수 설정 (MATLAB: ex1021.m 기준)
    T = 100                    # 초기 온도 (Initial Temperature)
    xl = np.array([-10, -10], dtype=float)  # 하한값 (Lower bounds)
    xu = np.array([10, 10], dtype=float)   # 상한값 (Upper bounds)
    rp = 0.8                   # 온도 감소 계수 (Reduction factor for T)
    rs = 0.9                   # 스텝 감소 계수 (Reduction factor for step size)
    crit = 1e-8                # 수렴 판정 기준 (Stopping criterion)

    try:
        # 3. 시뮬레이티드 어닐링 최적화 실행
        # xopt: 최적해, fopt: 최적 함수값, iter_count: 수렴 시점의 반복 횟수
        xopt, fopt, iter_count = saopt(fcy, T, xl, xu, rp, rs, crit)

        # 4. 결과 출력
        print("=== Simulated Annealing Optimization Result ===")
        print(f"최적해 (xopt): {xopt}")
        print(f"최적 함수값 (fopt): {fopt:.10f}")
        print(f"수렴 시 사이클 반복 횟수 (iter): {iter_count}")

    except Exception as e:
        print(f"최적화 실행 중 오류가 발생했습니다: {e}")

if __name__ == "__main__":
    main()