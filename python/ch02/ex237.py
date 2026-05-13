import numpy as np
from scipy.interpolate import RegularGridInterpolator

# 1. 데이터 정의
# T (온도): x축 데이터 (5개)
t_axis = np.array([300, 350, 400, 450, 500])
# P (압력): y축 데이터 (6개)
p_axis = np.array([150, 200, 225, 250, 275, 300])

# H (엔탈피) 데이터 테이블
# MATLAB의 H는 (P의 개수, T의 개수) 즉, (6, 5) 행렬입니다.
h_table = np.array([
    [3073.3, 3174.7, 3277.5, 3381.7, 3487.6],
    [3072.1, 3173.8, 3276.7, 3381.1, 3487.0],
    [3071.5, 3173.3, 3276.3, 3380.8, 3486.8],
    [3070.9, 3172.8, 3275.9, 3380.4, 3486.5],
    [3070.3, 3172.4, 3275.5, 3380.1, 3486.2],
    [3069.7, 3171.9, 3275.2, 3379.8, 3486.0]
])

# 2. 보간 함수 생성
# MATLAB의 'spline'과 대응하는 'cubic' 방법을 사용합니다.
# 입력 순서: (y축 좌표, x축 좌표) -> (P, T)
interp_h = RegularGridInterpolator((p_axis, t_axis), h_table, method='cubic')

# 3. 특정 지점 계산 (T=380, P=260)
# 주의: 함수 생성 시 (P, T) 순서로 정의했으므로 입력도 [260, 380] 순서입니다.
query_point = np.array([260, 380])
hi_result = interp_h(query_point)

print(f"Interpolated Enthalpy (Hi) at T=380, P=260: {hi_result[0]:.4f}")