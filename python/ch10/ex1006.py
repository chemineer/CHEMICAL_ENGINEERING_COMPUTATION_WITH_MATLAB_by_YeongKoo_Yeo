import numpy as np
from newtonopt import newtonopt

# 1. 목적 함수의 그레디언트(Gradient) 정의
# 뉴턴법(Newton's method)은 f(x) = 0인 지점을 찾으므로, 
# 최적화 문제에서는 목적 함수의 1차 도함수(그레디언트)를 전달해야 합니다.
def grad_fun(x):
    # 매트랩 원문: 100*(x(2)-x(1)^2)^2 + (1-x(1))^2
    # 위 식의 x[0]에 대한 편미분
    df_dx1 = -400 * x[0] * (x[1] - x[0]**2) - 2 * (1 - x[0])
    # 위 식의 x[1]에 대한 편미분
    df_dx2 = 200 * (x[1] - x[0]**2)
    return np.array([df_dx1, df_dx2])

# 2. 매개변수 설정
x0 = [-1.2, 1.0]   # 시작점
crit = 1e-6        # 허용 오차
kmax = 1000        # 최대 반복 횟수 (1e3)

# 3. newtonopt 함수 호출
# xopt: 최적점, fopt: 최적점에서의 그레디언트 값, iters: 수행된 반복 횟수
xopt, fopt, iters = newtonopt(grad_fun, x0, crit, kmax)

# 4. 결과 출력
print("--- Newton's Method 최적화 결과 ---")
print(f"최적점 (xopt): {xopt}")
print(f"최종 그레디언트 값 (fopt): {fopt}")
print(f"반복 횟수 (iter): {iters}")