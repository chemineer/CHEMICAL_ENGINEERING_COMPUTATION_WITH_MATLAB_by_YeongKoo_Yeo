# main.py 또는 실행 스크립트
from nnhzpipe import nnhzpipe

# 데이터 설정
L = 15
R = 0.009
K = 1e-6
n = 2
delP = 110

# 함수 호출
avgv = nnhzpipe(L, R, K, n, delP)

# 결과 출력
print(f"Average velocity (avgv) = {avgv[0]:.6f} m/s")