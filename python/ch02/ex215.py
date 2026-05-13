import numpy as np
from scipy.optimize import fsolve

# 1. 비선형 연립방정식 시스템 정의
def fun(x):
    # MATLAB: x(1), x(2), x(3) -> Python: x[0], x[1], x[2]
    f1 = np.sin(x[0]) + x[1]**2 + np.log(x[2]) - 7
    f2 = 3*x[0] + 2*x[1] - x[2]**3 + 1
    f3 = x[0] + x[1] + x[2] - 5
    return [f1, f2, f3]

# 2. 초기 추측값 설정 (x0 = [0, 2, 2])
x0 = [0, 2, 2]

# 3. fsolve를 사용하여 해 계산
solution = fsolve(fun, x0)

# 4. 결과 출력
print("--- 계산 결과 ---")
print(f"x1 = {solution[0]:.6f}")
print(f"x2 = {solution[1]:.6f}")
print(f"x3 = {solution[2]:.6f}")