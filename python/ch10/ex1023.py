import numpy as np
from ozopt import ozopt

def main():
    # 1. 문제 설정 (MATLAB: ex1023.m 기준)
    # nl: '<=' 제약 조건의 개수 (코드 내에서 '>='로 변환됨)
    # ne: '=' 등식 제약 조건의 개수
    nl = 2
    ne = 0
    
    # 제약 조건 행렬 A (각 행은 하나의 제약 조건을 나타내며, 마지막 열은 우변항 상수)
    # 2x1 - 5x2 <= 10  =>  -2x1 + 5x2 >= -10
    # 3x1 + 2x2 <= 9   =>  -3x1 - 2x2 >= -9
    A = np.array([
        [2, -5, 10],
        [3, 2, 9]
    ], dtype=float)
    
    # 목적 함수 계수 c (Minimize f = 4x1 + 3x2)
    c = np.array([4, 3], dtype=float)

    try:
        # 2. 0-1 정수 계획법 최적화 실행
        # xopt: 최적해 (0 또는 1의 값을 가짐)
        # fopt: 최적 목적 함수 값
        # iter_count: 반복 횟수
        xopt, fopt, iter_count = ozopt(A, c, nl, ne)

        # 3. 결과 출력
        print("=== Zero-One Programming Result ===")
        print(f"최적해 (xopt): {xopt}")
        print(f"최적 목적 함수 값 (fopt): {fopt}")
        print(f"반복 횟수 (iterations): {iter_count}")

    except Exception as e:
        print(f"최적화 실행 중 오류가 발생했습니다: {e}")

if __name__ == "__main__":
    main()