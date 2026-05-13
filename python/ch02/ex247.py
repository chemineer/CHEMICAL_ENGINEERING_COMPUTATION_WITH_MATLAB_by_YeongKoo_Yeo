import numpy as np
from eulerde import eulerde
from rk4th import rk4th

# 1. 미분 방정식 정의: dy/dt = f(t, y)
f = lambda t, y: 5 * np.exp(0.6 * t) - 2 * y

# 초기값 및 설정
y0 = 1.5
tspan = [0, 3]
n = 5
h = (tspan[1] - tspan[0]) / n  # 단계 크기(step size)

# --- 결과 출력 ---
print("--- Explicit Euler Method ---")
t_euler, y_euler = eulerde(f, tspan, y0, n)
print(np.hstack((t_euler, y_euler)))

print("\n--- 4th-order Runge-Kutta Method ---")
t_rk4, y_rk4 = rk4th(f, tspan, y0, n)
print(np.hstack((t_rk4, y_rk4)))