import numpy as np

def rk4th(f, tspan, y0, n):
    """
    4차 런지-쿠타 방법을 사용하여 dy/dt = f(t, y)를 해결합니다.
    
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
    
    # t 벡터 생성 (n+1개의 지점)
    t = np.linspace(t0, tf, n + 1)
    
    # y 벡터 초기화
    y = np.zeros(len(t))
    y[0] = y0
    
    # Runge-Kutta 4차 루프
    for k in range(n):
        # 현재 단계의 기울기들 계산
        k1 = f(t[k], y[k])
        k2 = f(t[k] + h/2, y[k] + h*k1/2)
        k3 = f(t[k] + h/2, y[k] + h*k2/2)
        k4 = f(t[k] + h, y[k] + h*k3)
        
        # 다음 단계 값 계산
        y[k+1] = y[k] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
        
    return t, y