import numpy as np
from rbopt import rbopt

def main():
    # 1. 목적 함수 정의 (MATLAB: f = @(x) x(1)^2 + 2*x(2)^2 + 2*x(1)*x(2))
    # 파이썬 인덱스는 0부터 시작하므로 x(1) -> x[0], x(2) -> x[1]로 매핑합니다.
    f = lambda x: x[0]**2 + 2 * x[1]**2 + 2 * x[0] * x[1]

    # 2. 초기값 및 설정값 (MATLAB: x0 = [0.5 1], crit = 1e-6)
    # x0: 시작점, crit: 수렴 판정 기준
    x0 = np.array([0.5, 1], dtype=float)
    crit = 1e-6

    try:
        # 3. Rosenbrock의 방법을 이용한 최적화 실행
        # xopt: 최적점, fopt: 최적점에서의 함수값, istag: 수행된 스테이지 수
        xopt, fopt, istag = rbopt(f, x0, crit)

        # 4. 결과 출력
        print("=== Rosenbrock's Method Optimization Result ===")
        print(f"최적해 (xopt): {xopt}")
        print(f"최적 함수값 (fopt): {fopt:.10f}")
        print(f"수행된 스테이지 수 (istag): {istag}")

    except Exception as e:
        print(f"최적화 실행 중 오류가 발생했습니다: {e}")

if __name__ == "__main__":
    main()