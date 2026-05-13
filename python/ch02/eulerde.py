import numpy as np

def eulerde(f, tspan, y0, n):
    """
    dy/dt = f(t, y)를 풀기 위한 명시적 오일러 방법 구현
    
    입력:
    f: 함수 dy/dt = f(t, y)
    tspan: [시작 시간, 종료 시간] 벡터
    y0: 초기값
    n: 구간의 개수
    
    출력:
    t: 독립 변수 벡터
    y: 종속 변수(해) 벡터
    """
    t0 = tspan[0]
    tf = tspan[1]
    h = (tf - t0) / n
    
    # t 벡터 생성 (t0부터 tf까지 n+1개의 지점)
    t = np.linspace(t0, tf, n + 1)
    
    # y 벡터 초기화
    y = np.zeros(len(t))
    y[0] = y0
    
    # 오일러 공식 적용: y(k) = y(k-1) + h * f(t(k-1), y(k-1))
    for k in range(1, len(t)):
        y[k] = y[k-1] + h * f(t[k-1], y[k-1])
        
    return t, y