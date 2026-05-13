import numpy as np
import matplotlib.pyplot as plt
from distBPeu import distBPeu  # distBPeu.py 모듈이 있다고 가정

# 1. 운영 데이터 및 파라미터 설정 (opdat)
opdat = {
    'N': 16,               # 총 단수
    'nc': 5,              # 성분 수
    'nf': np.array([5, 8]), # 피드 단 (Python 0-index 기준: MATLAB 6, 9 -> 5, 8)
    'fstate': 'V',        # 피드 상태
    'eos': 'pr',          # 상태 방정식 (Peng-Robinson)
    'criv': 1e-3          # 수렴 조건
}

# 배열 초기화
opdat['F'] = np.zeros(opdat['N'])
opdat['Tf'] = np.zeros(opdat['N'])
opdat['Pf'] = np.zeros(opdat['N'])
opdat['P'] = 240 * np.ones(opdat['N'])
opdat['V'] = np.zeros(opdat['N'])
opdat['L'] = np.zeros(opdat['N'])
opdat['U'] = np.zeros(opdat['N'])
opdat['W'] = np.zeros(opdat['N'])
opdat['Q'] = np.zeros(opdat['N'])
opdat['z'] = np.zeros((opdat['nc'], opdat['N']))

# 피드 및 초기 조건 설정
nf_idx = opdat['nf']
opdat['F'][nf_idx] = [41, 59]
opdat['V'][0] = 15
opdat['L'][0] = 150
opdat['U'][0] = 5
opdat['U'][2] = 3
opdat['W'][12] = 37  # MATLAB 13 -> Python 12

# 하부 유량 계산
opdat['L'][-1] = np.sum(opdat['F']) - opdat['V'][0] - opdat['U'][0] - opdat['U'][2] - opdat['W'][12]

opdat['Pf'][nf_idx] = [300, 275]
opdat['Tf'][nf_idx] = [170, 230]
opdat['Q'][2] = 2e5
opdat['z'][:, nf_idx[0]] = [0.061, 0.342, 0.463, 0.122, 0.012]
opdat['z'][:, nf_idx[1]] = [0.0085, 0.1017, 0.3051, 0.5085, 0.0762]

# 초기 온도 및 유량 추정값
opdat['T0'] = (300 / opdat['N']) * np.arange(1, opdat['N'] + 1)
opdat['V'][1:-1] = 170

# 2. 물성 데이터 설정 (mxdat)
mxdat = {
    'Pc': np.array([709.8, 617.4, 550.7, 489.5, 440.0]), # psia
    'Tc': np.array([550.0, 665.9, 765.3, 845.9, 914.2]), # R
    'k': np.zeros((opdat['nc'], opdat['nc'])),           # Binary interaction
    'w': np.array([0.1064, 0.1538, 0.1954, 0.2387, 0.2972]), # Acentric factors
    # Antoine parameters
    'Ant': np.array([
        [5.38389, 2847.921, 434.898],
        [5.35342, 3371.084, 414.488],
        [5.74162, 4126.385, 409.5179],
        [5.853654, 4598.287, 394.4148],
        [6.03924, 5085.758, 382.794]
    ]),
    # Specific heat parameters
    'Afi': np.array([
        [11.51606, 0.140309e-1, 0.0854034e-4, -0.110608e-7, 0.316220e-11],
        [15.58683, 0.2504953e-1, 0.1404258e-4, -0.352626e-7, 1.864467e-11],
        [20.79783, 0.314329e-1, 0.192851e-4, -0.458865e-7, 2.380972e-11],
        [25.64627, 0.389176e-1, 0.239729e-4, -0.584262e-7, 3.079918e-11],
        [30.17847, 0.519926e-1, 0.030488e-4, -0.27640e-7, 1.346731e-11]
    ])
}

# 3. BP Method 계산 호출
# distBPeu.py 내의 distBPeu 함수를 호출합니다.
x, y, T, L, V, iter_count = distBPeu(opdat, mxdat)

# 4. 결과 출력 및 시각화
print(f"Number of iterations (convergence criterion: {opdat['criv']}): {iter_count}")
print("Temperature (F):", T)

stages = np.arange(1, opdat['N'] + 1)

plt.figure(figsize=(12, 10))

# Temperature Plot
plt.subplot(2, 2, 1)
plt.plot(stages, T)
plt.xlabel('Stage')
plt.ylabel('T(F)')
plt.title('Temperature Profile')
plt.grid(True)

# Liquid Composition (x) Plot
plt.subplot(2, 2, 2)
markers = ['', '--', '-.', ':', '*']
for i in range(opdat['nc']):
    plt.plot(x[i, :], stages, label=f'x_{i+1}')
plt.xlabel('x')
plt.ylabel('Stage')
plt.gca().invert_yaxis() # 단수 위에서 아래로
plt.legend()
plt.title('Liquid Composition')
plt.grid(True)

# Vapor Composition (y) Plot
plt.subplot(2, 2, 3)
for i in range(opdat['nc']):
    plt.plot(y[i, :], stages, label=f'y_{i+1}')
plt.xlabel('y')
plt.ylabel('Stage')
plt.gca().invert_yaxis()
plt.legend()
plt.title('Vapor Composition')
plt.grid(True)

plt.tight_layout()
plt.show()