import numpy as np
from scipy.optimize import fsolve
from newtrapmv import *

# 1. 탈수소화 반응 시스템 함수 정의
def f(x):
    # MATLAB: K1 = 3.78, K2 = 0.137
    K1 = 3.78
    K2 = 0.137
    
    # x(1) -> x[0], x(2) -> x[1]
    # 분모와 분자가 복잡하므로 가독성을 위해 공통 항을 변수로 설정 가능합니다.
    x1, x2 = x[0], x[1]
    common_den1 = (1 - x1 - x2)
    common_den2 = (1 + x1 + 2*x2)
    numerator_term = (x1 + 2*x2)

    # 첫 번째 방정식: f1 = K1을 만족하는 식 (f1 - K1 = 0)
    f1 = (x1 * numerator_term) / (common_den1 * common_den2) - K1
    
    # 두 번째 방정식: f2 = K2를 만족하는 식 (f2 - K2 = 0)
    f2 = (x2 * (numerator_term**2)) / (common_den1 * (common_den2**2)) - K2
    
    return np.array([f1, f2])

# 2. 초기 추측값 및 실행
x0 = [0.80, 0.10]

# 이전 단계에서 정의한 newtrapmv 함수를 호출합니다.
try:
    solution = newtrapmv(f, x0)
    if solution is not None:
        print("\n--- 화학 평형 해석 결과 ---")
        print(f"반응 정도 x1: {solution[0]:.6f}")
        print(f"반응 정도 x2: {solution[1]:.6f}")
except NameError:
    print("오류: 'newtrapmv' 함수가 정의되지 않았습니다. 앞서 구현한 함수를 먼저 실행해주세요.")
    
# fsolve를 사용한 해결법
z_scipy = fsolve(f, x0)
print(f"SciPy fsolve 결과: {z_scipy}")