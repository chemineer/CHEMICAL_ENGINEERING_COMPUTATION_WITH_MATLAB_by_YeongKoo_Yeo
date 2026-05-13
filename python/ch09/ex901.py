import sympy

def calculate_laplace():
    # 1. 기호 변수 정의 (t: 시간 영역, s: 복소수 주파수 영역)
    # a, b는 상수로 정의
    t, s = sympy.symbols('t s')
    a, b = sympy.symbols('a b')

    # 2. 함수 f 정의
    # f = 1 + t + t^2 + sin(at) - t*cos(bt)
    f = 1 + t + t**2 + sympy.sin(a*t) - t*sympy.cos(b*t)

    # 3. 라플라스 변환 계산
    # sympy.laplace_transform(함수, 시간변수, 주파수변수)
    # 결과값은 (변환결과, 수렴조건, 수렴영역) 튜플로 반환되므로 [0]번 인덱스만 사용
    Lf = sympy.laplace_transform(f, t, s)[0]

    # 4. 결과 출력
    print("입력 함수 f(t):")
    sympy.pprint(f)
    print("\n라플라스 변환 결과 L{f(t)}:")
    sympy.pprint(Lf)

if __name__ == "__main__":
    calculate_laplace()