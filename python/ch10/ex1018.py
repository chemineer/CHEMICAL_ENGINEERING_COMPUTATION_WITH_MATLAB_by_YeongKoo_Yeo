import numpy as np
from hjopt import hjopt

def main():
    # 1. 목적 함수 정의 (MATLAB: f = @(x) 100*(x(1)^2 - x(2))^2 + (1 - x(1))^ 2)
    # 파이썬 인덱스는 0부터 시작하므로 x(1) -> x[0], x(2) -> x[1]로 매핑합니다.
    f = lambda x: 100 * (x[0]**2 - x[1])**2 + (1 - x[0])**2

    # 2. 초기값 및 설정값 (MATLAB: x0 = [-1 1], crit = 1e-6)
    x0 = np.array([-1, 1], dtype=float)
    crit = 1e-6

    try:
        # 3. Hooke-Jeeves 패턴 탐색 최적화 실행
        # xopt: 최적해, fopt: 최적 함수값, iter_count: 반복 횟수
        xopt, fopt, iter_count = hjopt(f, x0, crit)

        # 4. 결과 출력
        print("=== Hooke-Jeeves Pattern Search Result ===")
        print(f"최적해 (xopt): {xopt}")
        print(f"최적 함수값 (fopt): {fopt:.10f}")
        print(f"반복 횟수 (iterations): {iter_count}")

    except Exception as e:
        print(f"최적화 실행 중 오류가 발생했습니다: {e}")

if __name__ == "__main__":
    main()