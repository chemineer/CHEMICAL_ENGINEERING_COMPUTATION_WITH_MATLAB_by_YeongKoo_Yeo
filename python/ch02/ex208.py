from bisectn import bisectn
from secant import secant
from newtrap import newtrap

# 1. 함수 정의
# f(r) = 9496 * (1 - 12r^2 + 16r^3) - 1800
f = lambda r: 9.496e3 * (1 - 12*r**2 + 16*r**3) - 1800

# f'(r) = 9496 * (-24r + 48r^2)
df = lambda r: 9.496e3 * (-24*r + 48*r**2)

# 2. 수치해석 메서드 실행

# A. 이분법 (Bisection Method)
# 구간 [0, 0.5]에서 해를 찾음
sol_bisect = bisectn(f, 0, 0.5)
print(f"Bisection method: {sol_bisect:.6f}")

# B. 할선법 (Secant Method)
sol_secant = secant(f, 0, 0.5)
print(f"Secant method:    {sol_secant:.6f}")

# C. 뉴턴-랩슨법 (Newton-Raphson Method)
# 초기값 0.3과 도함수(fprime)를 함께 전달
sol_newton = newtrap(f, df, 0.3)
print(f"Newton-Raphson:   {sol_newton:.6f}")