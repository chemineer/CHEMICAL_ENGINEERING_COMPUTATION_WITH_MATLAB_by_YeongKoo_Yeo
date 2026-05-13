import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# 1. 데이터 정의
t = np.array([5, 10, 15, 20, 25, 30])  # 시간(hr)
x = np.array([0.245, 0.230, 0.211, 0.197, 0.187, 0.176])  # 자유 함수량
z = 0.03  # 나무 두께(m)

# 2. 비선형 모델 함수 정의
# C[0] = x0 (초기 수분 함량), C[1] = D (확산 계수)
def diffusion_model(t, x0, D):
    term_constant = (np.pi / (2 * z))**2
    # MATLAB 식: 8*x0/pi^2 * (exp(-D*t*term_constant) + exp(-9*D*t*term_constant)/9)
    term1 = np.exp(-D * t * term_constant)
    term2 = np.exp(-9 * D * t * term_constant) / 9
    return (8 * x0 / np.pi**2) * (term1 + term2)

# 3. 비선형 회귀 수행 (nlinfit 대응)
# 초기 추정값 C0 = [0.3, 0]
C0 = [0.3, 1e-6] # D는 보통 매우 작은 값이므로 0 대신 작은 양수를 권장합니다.
popt, pcov = curve_fit(diffusion_model, t, x, p0=C0)

x0_fit, D_fit = popt

# 4. 결과 출력
print(f"Initial free moisture content of the wood (x0) = {x0_fit:.6f}")
print(f"Diffusivity of water in the wood (D) = {D_fit:.6e}")

# 5. 시각화
tv = np.linspace(0, 30, 300)
xv = diffusion_model(tv, x0_fit, D_fit)

plt.figure(figsize=(8, 5))
plt.plot(tv, xv, label='Fitting by scipy curve_fit')
plt.plot(t, x, 'o', label='Data')
plt.grid(True)
plt.xlabel('t(hr)')
plt.ylabel('x(kg H2O / kg wood)')
plt.title('Estimation of Water Diffusivity in Wood')
plt.legend(loc='best')
plt.show()