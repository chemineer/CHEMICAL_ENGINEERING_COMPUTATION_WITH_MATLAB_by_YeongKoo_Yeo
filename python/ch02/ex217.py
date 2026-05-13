import numpy as np
from scipy import stats

# 1. 데이터 입력
tc = [
    5.394, 5.564, 5.654, 5.465, 5.495, 5.404, 5.524, 5.414,
    5.514, 5.335, 5.614, 5.455, 5.534, 5.524, 5.475, 5.295,
    5.384, 5.614, 5.554, 5.675, 5.455, 5.554, 5.504, 5.584
]

# 2. 통계량 계산
mean_val = np.mean(tc)
median_val = np.median(tc)
# 최빈값은 scipy.stats를 사용 (여러 개일 수 있으므로 첫 번째 값 선택)
mode_result = stats.mode(tc, keepdims=True)
mode_val = mode_result.mode[0]
# 분산과 표준편차 (MATLAB의 var, std와 동일하게 자유도 n-1 적용을 위해 ddof=1 설정)
var_val = np.var(tc, ddof=1)
std_val = np.std(tc, ddof=1)

# 3. 결과 출력
print(f"Mean               =  {mean_val:9.6f}")
print(f"Median             =  {median_val:9.6f}")
print(f"Mode               =  {mode_val:9.6f}")
print(f"Variance           =  {var_val:9.6f}")
print(f"Standard deviation =  {std_val:9.6f}")

import pandas as pd
df = pd.Series(tc)
print(df.describe())