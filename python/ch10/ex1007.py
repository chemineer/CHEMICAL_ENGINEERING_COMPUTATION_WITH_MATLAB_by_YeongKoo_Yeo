import numpy as np
from cgopt import cgopt

# 1. 목적 함수(fun) 및 그레디언트(delfun) 정의
# 매트랩: fun = @(x) 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
def fun(x):
    return 100 * (x[1] - x[0]**2)**2 + (1 - x[0])**2

# 매트랩: delfun = @(x) [-400*x(1)*(x(2)-x(1)^2)+2*(x(1)-1); 200*(x(2)-x(1)^2)];
def delfun(x):
    df_dx1 = -400 * x[0] * (x[1] - x[0]**2) + 2 * (x[0] - 1)
    df_dx2 = 200 * (x[1] - x[0]**2)
    return np.array([df_dx1, df_dx2])

# 2. 초기 매개변수 설정
x0 = [-1.2, 1.0]   # 시작점
alpha0 = 1.0       # 초기 스텝 사이즈
crit = 1e-6        # 중단 기준 (허용 오차)
kmax = 1000        # 최대 반복 횟수 (1e3)

# 3. 켤레구배법(cgopt) 함수 호출
xopt, fopt, iters = cgopt(fun, delfun, x0, alpha0, crit, kmax)

# 4. 결과 출력
print("--- 켤레구배법 (Conjugate Gradient) 최적화 결과 ---")
print(f"최적점 (xopt): {xopt}")
print(f"최솟값 (fopt): {fopt:.10f}")
print(f"반복 횟수 (iter): {iters}")