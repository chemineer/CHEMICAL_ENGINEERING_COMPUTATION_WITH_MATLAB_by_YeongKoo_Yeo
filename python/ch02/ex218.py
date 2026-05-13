import numpy as np
from scipy import stats
from scipy.integrate import quad

# 1. 데이터 설정
mu = 0       # 평균
sigma = 1    # 표준편차
lx = -0.8    # 하한값 (Lower bound)
ux = 0.8     # 상한값 (Upper bound)
zl = 2       # 기준 z-값

# 2. 확률 밀도 함수(PDF) 정의
# stats.norm.pdf(x, loc, scale)를 사용하여 정규 분포 함수 생성
fx = lambda x: stats.norm.pdf(x, mu, sigma)

# 3. 확률 계산

# A. 수치 적분(quad)을 이용한 방법
# MATLAB의 quad와 동일하게 PDF를 특정 구간에서 적분합니다.
P, error = quad(fx, lx, ux)

# B. 누적 분포 함수(CDF)를 이용한 방법
# Pr = P(lx <= Z <= ux) = CDF(ux) - CDF(lx)
Pr = stats.norm.cdf(ux, mu, sigma) - stats.norm.cdf(lx, mu, sigma)

# C. 특정 z-값보다 클 확률 (꼬리 영역)
# Pz = P(Z >= zl) = 1 - CDF(zl)
Pz = 1 - stats.norm.cdf(zl, mu, sigma)

# 4. 결과 출력
print(f'Probability of z assuming any value between {lx} and {ux} = {P:.6f}')
print(f'Probability of z lying between {lx} and {ux} = {Pr:.6f}')
print(f'Probability of observing z >= {zl} = {Pz:.6f}')