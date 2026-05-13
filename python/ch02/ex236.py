import numpy as np
from scipy.interpolate import RegularGridInterpolator

# 1. 데이터 정의
# RH (상대습도): x축 데이터
rh_axis = np.array([10, 30, 50, 70, 90]) 
# T (온도): 반드시 오름차순으로 변경 ([44, 51])
t_axis = np.array([44, 51])

# 데이터 테이블도 온도 축 순서에 맞춰 뒤집어줍니다. (44도 데이터가 먼저 오도록)
# 기존 h_table[0]이 51도, [1]이 44도였으므로 순서를 바꿉니다.
h_table = np.array([
    [6.52, 19.58, 32.70, 45.75, 58.78], # 44도 데이터
    [8.27, 24.75, 41.30, 58.10, 73.51]  # 51도 데이터
])

dp_table = np.array([
    [6.40, 23.34, 32.18, 38.32, 42.81], # 44도 데이터
    [10.10, 26.89, 37.01, 43.21, 47.80] # 51도 데이터
])

# 2. 보간 함수 생성 (RegularGridInterpolator)
# 데이터 포인트가 부족하므로 method를 'linear'로 수정합니다.
interp_h = RegularGridInterpolator((t_axis, rh_axis), h_table, method='linear')
interp_dp = RegularGridInterpolator((t_axis, rh_axis), dp_table, method='linear')

# 3. 특정 지점 계산 (RH=58.4, T=46.8)
# 입력 형식: [온도, 상대습도] 순서
point = np.array([46.8, 58.4])

hv_result = interp_h(point)
dpv_result = interp_dp(point)

# 결과 출력
print(f"Absolute Humidity (Hv): {hv_result[0]:.4f}")
print(f"Dew Point (DPv): {dpv_result[0]:.4f}")