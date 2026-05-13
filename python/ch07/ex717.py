# bindistMT.py 파일에서 함수를 불러옵니다.
from bindistMT import bindistMT

# 1. 데이터 입력 (Data)
alpha = 2.45    # 상대 휘발도
q = 0.80        # 원료의 열 상태 (q-line)
zf = 0.5        # 원료 농도
xd = 0.9        # 탑상 농도
xb = 0.1        # 탑저 농도
R = 1.5         # 환류비

# 2. 함수 호출
feedn, totaln = bindistMT(alpha, q, zf, xd, xb, R)

# 3. 결과 출력
print(f"Feed stage: {feedn}")
print(f"Total number of stages: {totaln}")