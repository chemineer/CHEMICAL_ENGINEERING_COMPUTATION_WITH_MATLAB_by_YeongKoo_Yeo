import numpy as np

# 1. 균일 분포(Uniform Distribution)에서 2x3 난수 생성
# MATLAB: x = rand(2,3)
x = np.random.rand(2, 3)
print("Uniform distribution (0 to 1):\n", x)

# 2. 표준 정규 분포(Normal Distribution)에서 2x3 난수 생성
# MATLAB: y = randn(2,3)
y = np.random.randn(2, 3)
print("\nNormal distribution (mean=0, std=1):\n", y)

# 3. 특정 범위 [3, 5] 사이의 3x3 난수 생성
# MATLAB: z = 3+(5-3)*rand(3,3)
z = 3 + (5 - 3) * np.random.rand(3, 3)
print("\nRandom numbers between 3 and 5:\n", z)