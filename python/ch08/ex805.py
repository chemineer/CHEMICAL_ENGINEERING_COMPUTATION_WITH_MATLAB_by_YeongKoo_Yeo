import numpy as np

# 1. 데이터 설정
T1 = 255
T4 = 298
L = np.array([0.015, 0.1, 0.075])    # 두께 (m)
k = np.array([0.151, 0.043, 0.762])  # 열전도도 (W/m.K)

# 2. 열유속(qx) 계산
# MATLAB의 L./k (요소별 나눗셈)는 numpy 배열의 / 연산과 동일합니다.
qx = (T1 - T4) / np.sum(L / k)

# 3. 결과 출력
print(f"qx = {qx:.4f}")