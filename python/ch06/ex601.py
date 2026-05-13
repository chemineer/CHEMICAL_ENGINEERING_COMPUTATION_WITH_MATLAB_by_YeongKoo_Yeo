import numpy as np
import matplotlib.pyplot as plt

# --- 데이터 설정 (source: 5) ---
T = np.array([313, 319, 323, 328, 333], dtype=float) # 온도
k = 1e-3 * np.array([0.43, 1.03, 1.80, 3.55, 7.17], dtype=float) # 속도 상수

# --- 아레니우스 식 분석 (ln(k) vs 1/T) ---
x = 1.0 / T # x축: 1/T
y = np.log(k) # y축: ln(k)

# polyfit(x, y, 1): 1차식(직선)으로 피팅[cite: 5]
# c[0]은 기울기(-E/R), c[1]은 절편(ln(A))
c = np.polyfit(x, y, 1)

# 빈도 인자(A)와 활성화 에너지(E) 계산
A = np.exp(c[1]) # A = exp(절편)[cite: 5]
E = -c[0] * 8.314 # E = -기울기 * R (R=8.314 J/mol·K)[cite: 5]

# 결과 출력
print(f"Pre-exponential factor (A): {A:g}")
print(f"Activation Energy (E): {E:g}")

# --- 시각화 ---
# 피팅된 직선을 그리기 위한 x값 생성[cite: 5]
xv = np.linspace(np.min(x), np.max(x), 100)
# polyval을 이용해 xv에 대응하는 y값(피팅값) 계산[cite: 5]
yv = np.polyval(c, xv)

plt.figure(figsize=(8, 6))
plt.plot(xv, yv, label='Fitted line') # 피팅 직선[cite: 5]
plt.plot(x, y, 'o', label='Data points') # 원본 데이터 포인트[cite: 5]

plt.xlabel('1/T(K^-1)') #[cite: 5]
plt.ylabel('ln(k)') #[cite: 5]
plt.title('Arrhenius Plot')
plt.legend()
plt.grid(True)
plt.show()