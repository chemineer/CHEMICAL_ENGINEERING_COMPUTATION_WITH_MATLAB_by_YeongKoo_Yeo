import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

def secpro(Kc=None):
    if Kc is None:
        Kc = float(input('Controller gain = '))
    
    # 시스템 정의: G(s) = 0.5*Kc / (0.5s^2 + s + 0.5*Kc)
    num = [0.5 * Kc]
    den = [0.5, 1, 0.5 * Kc]
    
    sys = signal.TransferFunction(num, den)
    t = np.linspace(0, 10, 101)
    t, y = signal.step(sys, T=t)
    
    plt.plot(t, y)
    plt.grid(True)
    plt.xlabel('Time t(sec)')
    plt.ylabel('Output C(t)')
    # plt.show()는 메인 루프가 끝난 후 호출하기 위해 제거함