# pipe_flow_calc.py
from hzpipe import hzpipe

# 1. 데이터 설정
L = 15.0       # 파이프 길이 (m)
R = 0.009      # 파이프 반지름 (m)
mu = 8.937e-4  # 점도 (Pa*s)
delP = 520.0   # 압력 강하 (Pa)

# 2. 함수 호출
Vavg = hzpipe(L, R, mu, delP)

# 3. 결과 출력
print(f"계산된 평균 유속 (Vavg): {Vavg:.6f} m/s")