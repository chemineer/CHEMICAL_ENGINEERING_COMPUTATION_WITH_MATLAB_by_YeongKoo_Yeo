import numpy as np
import math
from gaopt import gaopt

def main():
    # 1. 목적 함수 정의 (MATLAB: @fcy2)
    # f = cos(5*sqrt(sum(x.^2))) - 0.1*sum(x.^2) + 1
    # 유전 알고리즘은 이 값을 최대화하는 방향으로 탐색합니다.
    def fcy2(x):
        r2 = np.sum(x**2)
        term1 = math.cos(5 * math.sqrt(r2))
        term2 = 0.1 * r2
        return term1 - term2 + 1

    # 2. 매개변수 설정 (MATLAB: ex1022.m 기준)
    xl = np.array([-10, -10], dtype=float)  # 하한값 (Lower bounds)
    xu = np.array([10, 10], dtype=float)    # 상한값 (Upper bounds)
    nb = 8      # 변수당 비트 수 (Number of bits per variable)
    ps = 50     # 집단 크기 (Population size)
    ng = 60     # 세대 수 (Number of generations)
    mp = 0.05    # 변이 확률 (Mutation probability)

    try:
        # 3. 유전 알고리즘 최적화 실행
        # xopt: 최적해, fopt: 최적 함수값, total_evals: 총 함수 평가 횟수
        xopt, fopt, total_evals = gaopt(fcy2, xl, xu, nb, ps, ng, mp)

        # 4. 결과 출력
        print("=== Genetic Algorithm Optimization Result ===")
        print(f"최적해 (xopt): {xopt}")
        print(f"최적 함수값 (fopt): {fopt:.10f}")
        print(f"총 함수 평가 횟수: {total_evals}")

    except Exception as e:
        print(f"최적화 실행 중 오류가 발생했습니다: {e}")

if __name__ == "__main__":
    main()