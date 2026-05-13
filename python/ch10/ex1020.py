import numpy as np
from nmopt import nmopt

def main():
    # 1. 목적 함수 정의 (MATLAB: f = @(x) (x(1)+10*x(2))^2+5*(x(3)-x(4))^2+(x(2)-2*x(3))^4+10*(x(1)-x(4))^4)
    # 파이썬 인덱스는 0부터 시작하므로 x(1) -> x[0], x(2) -> x[1] 등으로 매핑합니다.
    f = lambda x: (x[0] + 10 * x[1])**2 + 5 * (x[2] - x[3])**2 + \
                  (x[1] - 2 * x[2])**4 + 10 * (x[0] - x[3])**4

    # 2. 초기값 및 설정값 (MATLAB: x0 = [-3 -1 0 1], crit = 1e-6)
    # x0: 초기 탐색 시작점, crit: 수렴 판정 기준
    x0 = np.array([-3, -1, 0, 1], dtype=float)
    crit = 1e-6

    try:
        # 3. Nelder-Mead 심플렉스 최적화 실행
        # xopt: 최적해, fopt: 최적 함수값, iter_count: 반복 횟수
        xopt, fopt, iter_count = nmopt(f, x0, crit)

        # 4. 결과 출력
        print("=== Nelder-Mead Simplex Optimization Result ===")
        print(f"최적해 (xopt): {xopt}")
        print(f"최적 함수값 (fopt): {fopt:.10f}")
        print(f"반복 횟수 (iterations): {iter_count}")

    except Exception as e:
        print(f"최적화 실행 중 오류가 발생했습니다: {e}")

if __name__ == "__main__":
    main()