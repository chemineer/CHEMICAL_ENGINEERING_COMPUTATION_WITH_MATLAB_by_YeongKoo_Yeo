import numpy as np
import matplotlib.pyplot as plt

# 1. Kc 범위 설정 (0부터 40까지 1씩 증가)
Kc_values = np.arange(0, 41, 1)

# 근을 저장할 리스트 초기화
all_roots = []

# 2. 각 Kc 값에 대한 다항식의 근 계산
for Kc in Kc_values:
    # 특성 방정식의 계수 정의: s^3 + 6s^2 + 11s + (6 + 2*Kc)
    coeff = [1, 6, 11, 6 + 2 * Kc]
    
    # 다항식의 근 계산 (MATLAB의 roots 함수 대응)
    sol = np.roots(coeff)
    all_roots.append(sol)

# 결과 데이터를 넘파이 배열로 변환
all_roots = np.array(all_roots)

# 3. 근궤적 시각화
plt.figure(figsize=(8, 6))

# 실수부(Real part)를 x축으로, 허수부(Imaginary part)를 y축으로 설정하여 '*' 모양으로 출력
plt.plot(all_roots.real, all_roots.imag, '*', markersize=5)

# 그래프 설정
plt.axhline(0, color='black', lw=1) # 실수축 강조
plt.axvline(0, color='black', lw=1) # 허수축 강조
plt.grid(True)
plt.xlabel('Real part')
plt.ylabel('Imaginary part')
plt.title('Root locus (Kc=0~40)')
plt.show()