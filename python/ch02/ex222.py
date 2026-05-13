import numpy as np
import matplotlib.pyplot as plt

# 1. 데이터 정의
t = np.array([1, 2, 3, 4, 5, 6, 7, 8])
T = np.array([50.8, 56.4, 55.1, 60.6, 61.5, 59.5, 54.1, 53.8])

# 그래프를 그리기 위한 정밀한 시간축 (1부터 8까지 0.1 간격)
tp = np.arange(1, 8.1, 0.1)

# 2. 다항식 적합 (Polynomial Fitting)
# polyfit(x, y, 차수) 순서로 입력합니다.
p1 = np.polyfit(t, T, 1) # 1차식
p2 = np.polyfit(t, T, 2) # 2차식
p4 = np.polyfit(t, T, 4) # 4차식

# 결과 계수 출력 (선택 사항)
print(f"1st-order coefficients: {p1}")
print(f"2nd-order coefficients: {p2}")
print(f"4th-order coefficients: {p4}")

# 3. 다항식 계산 (Calculate values)
T1 = np.polyval(p1, tp)
T2 = np.polyval(p2, tp)
T4 = np.polyval(p4, tp)

# 4. 결과 시각화
plt.figure(figsize=(10, 6))
plt.plot(t, T, 'o', label='Data')           # 원본 데이터
plt.plot(tp, T1, ':', label='1st-order')    # 1차 근사
plt.plot(tp, T2, '-.', label='2nd-order')   # 2차 근사 (.- 스타일 대응)
plt.plot(tp, T4, '-', label='4th-order')    # 4차 근사

plt.xlabel('t(hr)')
plt.ylabel('T(C)')
plt.grid(True)
plt.legend()
plt.title('Polynomial Fitting by polyfit function')
plt.show()