import numpy as np
import matplotlib.pyplot as plt

# 1. 실험 데이터 정의
S = np.array([1.2, 1.6, 3.2, 4.3, 5.8, 7.6, 8.8]) # 기질 농도 [S]
r = np.array([0.06, 0.12, 0.24, 0.27, 0.33, 0.34, 0.34]) # 반응 속도 r

# 2. Lineweaver-Burk Plot을 이용한 선형 회귀 (1/S와 1/r 관계)
# 1/r = (b/a)*(1/S) + (1/a) 형태의 1차식으로 적합
p = np.polyfit(1/S, 1/r, 1) # p[0]은 기울기(b/a), p[1]은 y절편(1/a)

# 3. 매개변수 추출
# a: 최대 반응 속도 (Vmax), b: 미하엘리스 상수 (Km)
a = 1 / p[1]
b = p[0] * a

print(f"추정된 매개변수 a (Vmax): {a:.4f}")
print(f"추정된 매개변수 b (Km): {b:.4f}")

# 4. 모델을 이용한 예측값 계산
Sv = np.arange(S[0], S[-1] + 0.1, 0.1) # S의 최소값부터 최대값까지 0.1 간격
rv = (a * Sv) / (b + Sv)

# 5. 시각화
plt.figure(figsize=(10, 6))
plt.plot(Sv, rv, label='Michaelis-Menten model') # 모델 예측 곡선
plt.plot(S, r, 'o', label='Experimental data')    # 실제 실험 데이터 점
plt.xlabel('[S]')
plt.ylabel('r')
plt.title('Michaelis-Menten Equation Fitting')
plt.legend()
plt.grid(True)
plt.show()