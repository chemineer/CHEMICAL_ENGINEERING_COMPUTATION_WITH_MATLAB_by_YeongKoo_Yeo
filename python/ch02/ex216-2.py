import numpy as np
from scipy.optimize import fsolve
from rxnfun import rxnfun  # rxnfun 함수가 정의된 rxnfun.py가 필요합니다.

def main():
    # 1. 초기 추정값 (x0) 및 온도(T) 설정
    x0 = np.array([0.001, 0.001, 0.001, 0.993, 1.0, 0.0001, 5.992, 1.0, 0.001, 10.0, 10.0, 10.0])
    T = 1000

    # 2. fsolve를 이용한 비선형 시스템 풀이
    # 매틀랩의 fsolve(@rxnfun, x0, [], T)와 같이 추가 파라미터 T를 args로 전달합니다.
    x = fsolve(rxnfun, x0, args=(T,))

    # 3. 라벨 정의
    comp = ['CH4', 'C2H4', 'C2H2', 'CO2', 'CO', 'O2', 'H2', 'H2O', 'C2H6']
    lamda = ['lambda1', 'lambda2', 'lambda3']

    # 4. 성분(Component) 결과 출력
    print('\n i\tComp. \t Initial Val.\t\tFinal val.')
    for k in range(len(comp)):
        # 매틀랩의 1-based 인덱스 출력을 위해 k+1 사용
        print(f'{k+1}\t{comp[k]:<6}\t{x0[k]:12.9f}\t{x[k]:12.9f}')

    # 5. 라그랑주 승수(Lambda) 결과 출력
    print('\n i\tLambda\tInitial Val.\tFinal val.')
    for i in range(len(lamda)):
        idx = len(comp) + i  # x 배열에서의 실제 인덱스
        print(f'{idx+1}\t{lamda[i]:<7}\t{x0[idx]:4.1f}\t{x[idx]:15.9f}')

if __name__ == "__main__":
    main()