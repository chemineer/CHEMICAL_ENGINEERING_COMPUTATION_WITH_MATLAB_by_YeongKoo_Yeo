import numpy as np
import matplotlib.pyplot as plt

def vtwall(rho, mu, delta):
    """
    수직 벽면을 따라 흐르는 액체 박막의 속도 분포 및 평균 속도 계산
    
    Parameters:
    rho   : 밀도 (kg/m^3)
    mu    : 점도 (kg/m/s)
    delta : 박막 두께 (m)
    
    Returns:
    avgv  : 평균 속도 (m/s)
    """
    g = 9.8
    # 100등분하여 위치 배열 생성
    h = delta / 100
    x = np.arange(0, delta + h, h)
    
    # 속도 분포 계산
    # Rg는 박막 유동의 특징적인 속도 스케일 (rho * g * delta^2 / 2 / mu)
    Rg = rho * g * delta**2 / (2 * mu)
    v = Rg * (1 - (x / delta)**2)
    
    # 평균 속도 계산
    avgv = Rg * 2 / 3
    
    # 그래프 시각화
    plt.figure()
    plt.plot(x, v)
    plt.grid(True)
    plt.ylabel('v_z(x)')
    plt.xlabel('x')
    plt.axis([0, delta, 0, np.max(v)])
    plt.show()
    
    return avgv