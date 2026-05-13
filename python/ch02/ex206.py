import numpy as np

# 1. 다항식의 계수 정의 (최고차항부터 상수항 순서)
# f(x) = 1x^5 - 3x^4 + 3x^3 - 2x^2 - 4x + 1
c = [1, -3, 3, -2, -4, 1]

# 2. 다항식의 근 계산
x = np.roots(c)

# 3. 결과 출력
print("The roots of the polynomial are:")
print(x)

# 가독성을 위해 근의 형태(실수/복소수)별로 출력해보기
for i, root in enumerate(x):
    if np.iscomplex(root):
        print(f"Root {i+1}: {root:.4f} (Complex)")
    else:
        print(f"Root {i+1}: {root.real:.4f} (Real)")
        
# 근을 대입하여 결과가 0에 수렴하는지 확인
check = np.polyval(c, x)
print("\nVerification (should be close to 0):")
print(check)