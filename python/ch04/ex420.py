import numpy as np
from scipy.optimize import fsolve

# --- 사용자 정의 함수들 ---

def gamma12(x1, T):
    """활량 계수 계산: x1이 배열인 경우를 대비해 처리"""
    # x1이 배열이면 첫 번째 요소만 추출하여 계산
    if isinstance(x1, (np.ndarray, list)):
        x1 = x1[0]
        
    A_val = 2.771 - 0.00523 * T
    gamma1 = np.exp(A_val * (1 - x1)**2)
    gamma2 = np.exp(A_val * x1**2)
    return np.array([gamma1, gamma2])

def vp12(A, B, C, T):
    """Vapor pressure by Antoine equation"""
    Psat = np.exp(A - B / (T - C))
    return Psat

def dewpf(x1, y, T, Psat):
    """(2) Dew P calculation"""
    # fsolve에서 넘어온 x1이 배열일 경우를 대비해 float 추출
    x1_val = float(x1[0]) if hasattr(x1, "__len__") else float(x1)
    x = np.array([x1_val, 1 - x1_val])
    gamma = gamma12(x1_val, T)
    P = np.sum(x * gamma * Psat)
    # 결과를 스칼라로 반환하기 위해 .item() 사용
    res = np.sum(y * P / (gamma * Psat)) - 1
    return res.item() if hasattr(res, "item") else res

def bubtf(T, x1, A, B, C, P):
    """(3) Bubble T calculation"""
    T_val = float(T[0]) if hasattr(T, "__len__") else float(T)
    x = np.array([x1, 1 - x1])
    gamma = gamma12(x1, T_val)
    Psat = vp12(A, B, C, T_val)
    res = np.sum(x * gamma * Psat) - P
    return res.item() if hasattr(res, "item") else res

def dewfun(T, y1, x1, A, B, C, P):
    """(4) Dew T calculation"""
    T_val = float(T[0]) if hasattr(T, "__len__") else float(T)
    y = np.array([y1, 1 - y1])
    Psat = vp12(A, B, C, T_val)
    gamma = gamma12(x1, T_val)
    res = np.sum(y * P / (gamma * Psat)) - 1
    return res.item() if hasattr(res, "item") else res

# --- 실행 부분 ---

A_ant = np.array([16.59158, 14.25326])
B_ant = np.array([3643.31, 2665.54])
C_ant = np.array([33.424, 53.424])

# (1) Bubble P
print("--- (1) Bubble P ---")
T, x1_val = 318.15, 0.25
x = np.array([x1_val, 1 - x1_val])
gamma = gamma12(x1_val, T)
Psat = vp12(A_ant, B_ant, C_ant, T)
P = np.sum(x * gamma * Psat)
y = (x * gamma * Psat) / P
print(f"P = {P:.4f}, y1 = {y[0]:.4f}, y2 = {y[1]:.4f}\n")
print(f"gamma1 = {gamma[0]:.5f}, gamma2 = {gamma[1]:.5f}\n")

# (2) Dew P
print("--- (2) Dew P ---")
T, y1_val = 318.15, 0.60
y = np.array([y1_val, 1 - y1_val])
Psat = vp12(A_ant, B_ant, C_ant, T)
x1_sol = fsolve(dewpf, 0.7, args=(y, T, Psat))[0]
x = np.array([x1_sol, 1 - x1_sol])
gamma = gamma12(x1_sol, T)
P = np.sum(x * gamma * Psat)
print(f"P = {P:.4f}, x1 = {x[0]:.4f}, x2 = {x[1]:.4f}\n")
print(f"gamma1 = {gamma[0]:.5f}, gamma2 = {gamma[1]:.5f}\n")

# (3) Bubble T
print("--- (3) Bubble T ---")
P_target, x1_val = 101.33, 0.85
x = np.array([x1_val, 1 - x1_val])
T_sol = fsolve(bubtf, 300.0, args=(x1_val, A_ant, B_ant, C_ant, P_target))[0]
gamma = gamma12(x1_val, T_sol)
Psat = vp12(A_ant, B_ant, C_ant, T_sol)
y_sol = (x * gamma * Psat) / P_target
print(f"T = {T_sol:.4f}, y1 = {y_sol[0]:.4f}, y2 = {y_sol[1]:.4f}\n")
print(f"gamma1 = {gamma[0]:.5f}, gamma2 = {gamma[1]:.5f}\n")

# (4) Dew T
print("--- (4) Dew T ---")
P_target, y1_val = 101.33, 0.40
y = np.array([y1_val, 1 - y1_val])
T_guess, x1_guess = 330.0, 0.5 # 초기값
crx, cx = 1.0, 1e-6

while crx > cx:
    # fsolve 결과의 첫 번째 요소[0]를 가져와 스칼라로 사용
    T_guess = fsolve(dewfun, T_guess, args=(y1_val, x1_guess, A_ant, B_ant, C_ant, P_target))[0]
    gamma = gamma12(x1_guess, T_guess)
    Psat = vp12(A_ant, B_ant, C_ant, T_guess)
    x_new = y * P_target / (gamma * Psat)
    crx = np.abs(x_new[0] - x1_guess)
    x1_guess = x_new[0]
    x_final = x_new

print(f"T = {T_guess:.4f}, x1 = {x_final[0]:.4f}, x2 = {x_final[1]:.4f}")
print(f"gamma1 = {gamma[0]:.5f}, gamma2 = {gamma[1]:.5f}\n")