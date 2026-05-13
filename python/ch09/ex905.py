import control as ctrl
import numpy as np

# 분자(num)와 분모(den) 다항식 정의
# num = 2s + 1
num = [2, 1]

# 분모: (s+1)(s+2) = s^2 + 3s + 2
# conv 함수는 다항식 곱셈을 수행합니다.
den = np.convolve([1, 1], [1, 2])

# 전달함수 G(s) 생성
G = ctrl.tf(num, den)

# 결과 출력
print("전달함수 G(s):")
print(G)