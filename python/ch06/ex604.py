import numpy as np
from scipy.optimize import fsolve

# --- 평형 방정식 정의 ---
# MATLAB: f = @(x) 148.4 - x^2/(1 - x)^2
def equilibrium_func(x):
    # 평형 상수 K = [CO2][H2] / ([CO][H2O]) 식을 정리한 형태
    return 148.4 - (x**2) / (1 - x)**2

# --- 수치 해법 실행 ---
# 초기 추정값 x0 = 0.5
x0 = 0.5

# fsolve를 사용하여 방정식의 해(root)를 찾습니다.
# MATLAB의 fzero와 대응되는 라이브러리입니다.
solution = fsolve(equilibrium_func, x0)

# fsolve는 배열을 반환하므로 첫 번째 요소를 추출합니다.
x_val = solution[0]

# --- 결과 출력 ---
print(f"Equilibrium conversion (x): {x_val:.6f}")