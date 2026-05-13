import numpy as np
from scipy.optimize import fsolve
from newtrapmv import newtrapmv

# 1. 비선형 연립방정식 시스템 정의
def f(x):
    # x(1) -> x[0], x(2) -> x[1], x(3) -> x[2]
    f1 = np.cos(x[0]) + x[1]**2 + np.log(x[2]) - 8
    f2 = 4*x[0] + 3**x[1] - x[2]**3 + 2
    f3 = x[0] + x[1] + x[2] - 6
    return np.array([f1, f2, f3])

# 2. 초기 추측값 설정
x0 = np.array([1, 1, 1])

# 3. 해 구하기 (앞서 작성한 newtrapmv 함수 활용)
# 만약 함수가 같은 파일에 있다면 바로 호출 가능합니다.
try:
    z = newtrapmv(f, x0)
    print("--- 최종 결과 ---")
    print(f"x1: {z[0]:.6f}")
    print(f"x2: {z[1]:.6f}")
    print(f"x3: {z[2]:.6f}")
except Exception as e:
    print(f"오류 발생: {e}")
    

# fsolve를 사용한 해결법
z_scipy = fsolve(f, x0)
print(f"SciPy fsolve 결과: {z_scipy}")