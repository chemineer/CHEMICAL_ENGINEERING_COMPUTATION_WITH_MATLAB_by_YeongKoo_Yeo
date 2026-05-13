import control as ct
import matplotlib.pyplot as plt
import numpy as np

# 1. 주파수 범위 설정 (10^-2 ~ 10^1 rad/min, 300개 지점)
# 매트랩: w = logspace(-2, 1, 300);
w = np.logspace(-2, 1, 300)

# 2. 감쇠비(zeta) 설정 (0부터 1까지 0.25 간격)
# 매트랩: zeta = [0:0.25:1];
zeta_list = np.arange(0, 1.25, 0.25)

num = [1]

plt.figure(figsize=(10, 8))

for i, zeta in enumerate(zeta_list):
    den = [2.25, 3 * zeta, 1]
    sys = ct.tf(num, den)
    
    # frequency_response는 대부분의 버전에서 안정적으로 작동합니다.
    response = ct.frequency_response(sys, omega=w)
    
    # 데이터 추출
    mag = response.magnitude
    phase = response.phase  # 이 값은 라디안 단위입니다.
    
    # 5. 그래프 그리기
    # Magnitude Plot (상단)
    plt.subplot(2, 1, 1)
    plt.loglog(w, mag, label=f'zeta={zeta}')
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.ylabel('AR (Amplitude Ratio)')
    plt.title('Response of 2nd-order process')
    
    # Phase Plot (하단)
    plt.subplot(2, 1, 2)
    plt.semilogx(w, np.degrees(phase)) # 도(degree) 단위로 변환하여 출력
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.xlabel('w (rad/min)')
    plt.ylabel('Phase (deg)')

# 6. 매트랩의 text 함수와 유사하게 특정 위치에 설명 추가
plt.subplot(2, 1, 1)
plt.text(0.8, 20, 'zeta=0')
plt.text(0.4, 0.2, 'zeta=1')
plt.legend(loc='best')

plt.tight_layout()
plt.show()