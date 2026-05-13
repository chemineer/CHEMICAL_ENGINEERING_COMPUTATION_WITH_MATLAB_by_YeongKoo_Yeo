import control as ct
import matplotlib.pyplot as plt
import numpy as np

# 1. 시스템 정의: G(s) = 0.8 / (5s + 1)
num = [0.8]
den = [5, 1]
delay = 2

# 2. 방법 1: 주파수 응답 데이터를 직접 수정하여 마진 계산
# Bode 데이터를 추출한 뒤 시간 지연에 따른 위상 지연을 수동으로 반영
G_base = ct.tf(num, den)
# 주파수 범위 설정 (충분한 범위를 확보)
w = np.logspace(-2, 1, 1000)
mag, phase, omega = ct.bode(G_base, w, plot=False)

# 시간 지연에 의한 위상 변화 반영: Phase_total = Phase_base - (delay * omega)
# ct.bode의 phase는 라디안 단위이므로 그대로 계산 후 필요시 변환
phase_delayed = phase - (delay * omega)

# 이득 여유(Gm), 위상 여유(Pm) 등 계산
# ct.margin은 데이터(mag, phase, omega)를 직접 입력받을 수 있습니다.
Gm, Pm, Wcg, Wcp = ct.margin(mag, np.degrees(phase_delayed), omega)

print(f"Manual Calculation Results:")
print(f"Gain Margin: {Gm}, Phase Margin: {Pm}")
print(f"Wcg: {Wcg}, Wcp: {Wcp}\n")

# 3. 방법 2: 시스템 객체에 지연 시간을 설정하여 시각화 (권장 방식)
G_delayed = ct.tf(num, den)
G_delayed.ioTimeDelay = delay

plt.figure(figsize=(10, 8))
# ct.stability_margins 또는 ct.bode_plot을 사용하면 마진이 표시된 그래프를 얻을 수 있습니다.
# 최신 버전에서는 ct.margin_plot을 사용하여 매트랩의 margin(G)과 유사한 결과를 냅니다.
ct.bode_plot(G_delayed, omega=w, margins=True)

plt.suptitle('Bode Plot with Stability Margins (Example 9.28)')
plt.show()