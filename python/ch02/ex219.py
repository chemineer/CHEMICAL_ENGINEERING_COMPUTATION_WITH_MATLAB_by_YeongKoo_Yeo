import numpy as np
from scipy import stats

# 1. 데이터 입력
abdat = np.array([0.72, 0.54, 0.62, 0.80, 0.76, 0.64, 0.75, 0.94, 0.85, 0.44])

# 2. 기초 통계량 계산
avg = np.mean(abdat)
stdv = np.std(abdat, ddof=1) # MATLAB의 std와 맞추기 위해 자유도 1 설정

# 3. t-검정 수행 (귀무가설 mu = 0.86)
mu = 0.86
alpha = 0.05 # 유의수준

# ttest_1samp 함수는 기본적으로 양측 검정(both tails)을 수행합니다
t_stat, p_value = stats.ttest_1samp(abdat, mu)

# 4. 가설 채택/기각 여부 결정 (h)
# p-value가 유의수준(alpha)보다 작으면 귀무가설 기각 (h=1)
h = 1 if p_value < alpha else 0

# 5. 결과 출력
print(f'h (result of the hypothesis test) = {h}')
print(f'p (probability) = {p_value:.6f}')