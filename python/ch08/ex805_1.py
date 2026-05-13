import numpy as np

# 1. 데이터 설정
T1 = 255
T4 = 298
L = np.array([0.015, 0.075])  # 기존 층들의 두께
k = np.array([0.151, 0.762])  # 기존 층들의 열전도도
kB = 0.043                    # 추가할 재료(B)의 열전도도

# 2. dLB (재료 B의 두께) 계산
# MATLAB의 L./k는 numpy의 / 연산과 동일하며, sum()은 np.sum()으로 대응됩니다.
dLB = kB * ((T1 - T4) / (-8.5204) - np.sum(L / k))

# 3. 결과 출력
print(f"dLB = {dLB:.4f}")