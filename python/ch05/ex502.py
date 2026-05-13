# annulus_flow_calc.py
from hzannpipe import hzannpipe

# 1. 데이터 설정
mu = 8.937e-4   # 점도 (Pa*s)
delP = 110.0    # 압력 강하 (Pa)
L = 12.0        # 관 길이 (m)
R1 = 0.025      # 내부 반경 (m)
R2 = 0.036      # 외부 반경 (m)

# 2. 함수 호출
avgv = hzannpipe(L, R1, R2, mu, delP)

# 3. 결과 출력
print(f"환형 관 내 평균 유속 (avgv): {avgv:.6f} m/s")