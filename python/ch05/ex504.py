# terminal_velocity.py
from scipy.optimize import fsolve
from vtfun import vtfun

# 1. 데이터 설정
rp = 1780.0      # 입자 밀도 (kg/m^3)
ro = 994.6       # 유체 밀도 (kg/m^3)
dp = 2e-4        # 입자 직경 (m)
mu = 8.931e-4    # 유체 점도 (Pa*s)
vt0 = 1e-3       # 초기 추정치 (m/s)

# 2. fsolve를 이용한 방정식 풀이
# MATLAB의 fzero(@vtfun, vt0, [], rp, ro, mu, dp)에 대응
vt = fsolve(vtfun, vt0, args=(rp, ro, mu, dp))

print(f"입자의 종말 속도 (vt): {vt[0]:.6f} m/s")