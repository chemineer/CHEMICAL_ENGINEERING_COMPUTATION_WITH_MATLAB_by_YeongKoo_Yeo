import numpy as np
from scipy.integrate import trapezoid, cumulative_trapezoid
#numpy 1.25 버전부터 np.trapz가 Deprecated(사용 중단 예정) 되었고, 최신 버전인 numpy 2.0부터는 아예 삭제됨

# 1. 데이터 정의
t = np.array([0, 0.5, 1.2, 1.6, 2.5, 3.1, 4.8, 6.9])
v = np.array([0, 5, 12, 15, 23, 28, 38, 47])

t1 = np.array([0, 0.5, 1.2, 1.6, 2.5])
v1 = np.array([0, 5, 12, 15, 23])

# 2. trapezoid: 전체 면적 계산 (기존 np.trapz 대체)
# scipy.integrate.trapezoid(y, x) 순서로 입력합니다.
z1 = trapezoid(v1, t1)

# 3. cumulative_trapezoid: 누적 면적 계산 (기존 cumtrapz와 동일)
# initial=0을 설정해야 MATLAB의 cumtrapz와 결과 배열 길익가 같아집니다.
z = cumulative_trapezoid(v, t, initial=0)

# 4. 결과 출력
print(f"z1 결과 (전체 적분): {z1:.2f}")
print("-" * 40)
print("z (누적 적분) 결과:")
print(z)