import numpy as np
from bnbopt import bnbopt

def main():
    # 1. 문제 설정 (MATLAB: ex1025.m 기준)
    
    # 제약 조건 계수 행렬 A (6x9 행렬)
    A = np.array([
        [1, 0, 0, 1, 0, 1, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0, 1],
        [1, 1, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 1]
    ], dtype=float)
    
    # 우변항 상수 벡터 b (모두 1인 1x6 벡터)
    b = np.ones(6, dtype=float)
    
    # 목적 함수 계수 벡터 c (Minimize f = -[5 3 1 3 5 2 5 5 2])
    # MATLAB: c = -[5 3 1 3 5 2 5 5 2]
    c = -np.array([5, 3, 1, 3, 5, 2, 5, 5, 2], dtype=float)
    
    # 제약 조건의 종류별 개수
    nl = 3  # '<=' 제약 조건의 개수 (처음 3개)
    ne = 3  # '=' 등식 제약 조건의 개수 (나중 3개)
    ng = 0  # '>=' 제약 조건의 개수
    
    # 정수 제약이 있는 변수의 인덱스 (1-based index)
    # ibd = [1:9] -> 모든 변수(x1 ~ x9)가 정수여야 함
    ibd = np.arange(1, 10) 

    try:
        # 2. 분기한정법 최적화 실행
        # xopt: 최적해, fopt: 최적 목적 함수 값, iter_count: 반복 횟수
        xopt, fopt, iter_count = bnbopt(A, b, c, nl, ng, ne, ibd)

        # 3. 결과 출력
        print("=== Branch and Bound Optimization Result (Ex 10.25) ===")
        if xopt is not None:
            print(f"최적해 (xopt):\n{xopt}")
            print(f"최적 목적 함수 값 (fopt): {fopt:.4f}")
            print(f"반복 횟수 (iterations): {iter_count}")
        else:
            print("최적해를 찾지 못했습니다.")

    except Exception as e:
        print(f"최적화 실행 중 오류가 발생했습니다: {e}")

if __name__ == "__main__":
    main()