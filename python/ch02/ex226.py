import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 1. 데이터 정의
t = np.array([2.1, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0, 3.1])
v = np.array([3.091, 2.699, 1.801, 1.698, 1.412, 1.297, 0.702, 0.597])
z = np.linspace(min(t), max(t), 100)

# --- 방법 1: 선형화 모델을 이용한 Least Squares (행렬 연산) ---
# 모델: v = A * exp(-t) / (1 + B*t) => 1/v = (1/A)*exp(t) + (B/A)*t*exp(t)
X = np.column_stack([np.exp(t), t * np.exp(t)])
Y = 1 / v
C = np.linalg.inv(X.T @ X) @ X.T @ Y

A_LS = 1 / C[0]
B_LS = A_LS * C[1]
yz_LS = A_LS * np.exp(-z) / (1 + B_LS * z)

print(f"방법 1 (Least Squares) 결과: A = {A_LS:.4f}, B = {B_LS:.4f}")

# --- 방법 2: Scipy의 curve_fit (nlinfit 대응) ---
# 사용할 비선형 모델 함수 정의
def model_func(t, A, B):
    return A * np.exp(-t) / (1 + B * t)

# 초기값 (D0 = [1, 1]) 설정 및 최적화
D0 = [1, 1]
popt, _ = curve_fit(model_func, t, v, p0=D0)

A_fit, B_fit = popt
yz_fit = model_func(z, A_fit, B_fit)

print(f"방법 2 (curve_fit) 결과: A = {A_fit:.4f}, B = {B_fit:.4f}")

# 3. 시각화
plt.figure(figsize=(9, 6))
plt.plot(z, yz_LS, label='Nonlinear regression (LS)')
plt.plot(z, yz_fit, '--', label='Built-in fun curve_fit')
plt.plot(t, v, 'o', label='Data')

plt.xlabel('t')
plt.ylabel('v')
plt.title('Nonlinear Regression Comparison')
plt.legend(loc='best')
plt.grid(True)
plt.show()