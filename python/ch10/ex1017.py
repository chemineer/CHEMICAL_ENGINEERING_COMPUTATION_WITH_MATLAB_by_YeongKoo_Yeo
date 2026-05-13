import numpy as np
from cycopt import cycopt

def main():
    # 1. 목적 함수 정의 (MATLAB의 fc 수식을 파이썬 함수로 변환)
    # fc = @(x) (x(1)+10*x(2))^2+5*(x(3)-x(4))^2+(x(2)-2*x(3))^4+10*(x(1)-x(4))^4
    def fc(x):
        # 파이썬 인덱스는 0부터 시작함에 유의
        term1 = (x[0] + 10 * x[1])**2
        term2 = 5 * (x[2] - x[3])**2
        term3 = (x[1] - 2 * x[2])**4
        term4 = 10 * (x[0] - x[3])**4
        return term1 + term2 + term3 + term4

    # 2. 초기값 및 설정값 (MATLAB: x0 = [-3 -1 0 1], crit = 1e-4)
    x0 = np.array([-3, -1, 0, 1], dtype=float)
    crit = 1e-4

    try:
        # 3. 순차적 좌표 탐색(Cyclic Coordinate Search) 실행
        # xopt: 최적해, fopt: 최적 함수값, iter_count: 반복 횟수
        xopt, fopt, iter_count = cycopt(fc, x0, crit)

        # 4. 결과 출력
        print("=== Cyclic Coordinate Search Result ===")
        print(f"최적해 (xopt): {xopt}")
        print(f"최적 함수값 (fopt): {fopt:.10f}")
        print(f"반복 횟수 (iterations): {iter_count}")

    except Exception as e:
        print(f"최적화 실행 중 오류가 발생했습니다: {e}")

if __name__ == "__main__":
    main()