import sympy as sp

# 심볼 정의
s = sp.symbols('s')
t = sp.symbols('t', real=True, positive=True)

# 함수 F(s) 정의
F = 2*s / (s**2 + 4*s + 1)

# 라플라스 역변환 수행 (F를 s에서 t로 변환)
f = sp.inverse_laplace_transform(F, s, t)

# 결과 출력
print(f.simplify())