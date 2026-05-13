import numpy as np
import matplotlib.pyplot as plt

# --- 데이터 설정 ---
t = np.array([0, 2, 4, 7, 9, 18], dtype=float)
Ca = np.array([1.48, 1.01, 0.67, 0.58, 0.51, 0.32], dtype=float)

# --- 행렬 연산을 통한 선형 회귀 (Least Squares) ---
# M = [t' ones]
M = np.column_stack([t, np.ones(len(t))])
# Mi = inv(M'*M)*M' (의사 역행렬 계산)
Mi = np.linalg.inv(M.T @ M) @ M.T

# 각 차수별 종속 변수 설정
Y0 = Ca           # 0차: Ca = -kt + Ca0
Y1 = np.log(Ca)    # 1차: ln(Ca) = -kt + ln(Ca0)
Y2 = 1.0 / Ca     # 2차: 1/Ca = kt + 1/Ca0

# 계수 계산 (X = Mi * Y)
X0 = Mi @ Y0
X1 = Mi @ Y1
X2 = Mi @ Y2

# 속도 상수(k) 및 초기 농도(Ca0) 추출[cite: 7]
k0 = -X0[0]
Ca0_res = X0[1]

k1 = -X1[0]
Ca1_res = np.exp(X1[1])

k2 = X2[0]
Ca2_res = 1.0 / X2[1]

# 결과 출력[cite: 7]
print(f"0th order: k = {k0:g}, Ca0 = {Ca0_res:g}")
print(f"1st order: k = {k1:g}, Ca0 = {Ca1_res:g}")
print(f"2nd order: k = {k2:g}, Ca0 = {Ca2_res:g}")

# --- 시각화를 위한 계산 ---[cite: 7]
tv = np.arange(0, 20.1, 0.1) # 0부터 20까지 0.1 간격[cite: 7]

Caz = -k0 * tv + Ca0_res                      # 0차 예측값[cite: 7]
Caf = np.exp(-k1 * tv + np.log(Ca1_res))      # 1차 예측값[cite: 7]
Cas = 1.0 / (k2 * tv + 1.0 / Ca2_res)         # 2차 예측값[cite: 7]

# --- 그래프 그리기 ---[cite: 7]
plt.figure(figsize=(10, 6))
plt.plot(tv, Caz, ':', label='0th order')     # 점선[cite: 7]
plt.plot(tv, Caf, '--', label='1st order')    # 파선[cite: 7]
plt.plot(tv, Cas, '-', label='2nd order')     # 실선[cite: 7]
plt.plot(t, Ca, 'o', label='Data')            # 원본 데이터[cite: 7]

plt.xlabel('t(sec)')
plt.ylabel('C_A(mol/l)')
plt.legend(loc='best')
plt.title('Determination of Reaction Order')
plt.grid(True)
plt.show()