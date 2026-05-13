# falling_film_calc.py
from vtwall import vtwall

# 1. 데이터 설정
rho = 780.0    # 밀도 (kg/m^3)
mu = 0.172     # 점도 (Pa*s)
delta = 0.0021 # 액막 두께 (m)

# 2. 함수 호출
avgv = vtwall(rho, mu, delta)

# 3. 결과 출력
print(f"계산된 액막의 평균 유속 (avgv): {avgv:.6f} m/s")