# calculate_friction.py
from frfactor import frfactor

# 1. 데이터 설정 및 호출
# 첫 번째 케이스
eD1, Nre1 = 0, 2e4
frfactor(eD1, Nre1)

# 두 번째 케이스
eD2, Nre2 = 0, 3.2e7
frfactor(eD2, Nre2)