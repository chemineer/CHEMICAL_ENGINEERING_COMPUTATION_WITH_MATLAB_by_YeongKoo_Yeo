import numpy as np

# 1. 상수 및 범위 설정
D = 3.26e-8  # 확산 계수
# x 범위: -1e-3부터 1e-3까지 1e-4 간격
x = np.arange(-1e-3, 1.1e-3, 1e-4) 

# 2. 농도 함수 C(x) 정의
C = lambda x: -1.3e6*x**4 + 7.1e3*x**3 - 14*x**2 - 0.364*x + 0.001

# 3. 수치 미분 수행
f = C(x)
# diff(f) / diff(x)를 통해 농도 기울기(dC/dx) 계산
numf = np.diff(f) / np.diff(x)

# 4. CO2 플럭스 계산 (Fick의 확산 법칙: J = -D * dC/dx)
NCO2 = -D * numf

# 5. 양 끝단의 플럭스 및 순 플럭스(Net Flux) 추출
nfl = NCO2[0]       # 왼쪽 끝 (x = -0.001)
nfr = NCO2[-1]      # 오른쪽 끝 (x = 0.001)
netf = nfl - nfr    # 순 플럭스

# 6. 결과 출력
print(f"Flux at x = -0.001: {nfl:g}, Flux at x = 0.001: {nfr:g}, Net flux: {netf:g}")