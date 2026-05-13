import numpy as np
from compID import compID

def gfree(T, cname):
    """
    저압 상태 기체의 깁스 생성 자유 에너지(kcal/gmol) 추정 함수
    
    입력:
    - T: 온도 (K) (스칼라 또는 배열 가능)
    - cname: 화합물 이름 또는 화학식
    
    출력:
    - gf: 깁스 생성 자유 에너지 (kcal/gmol)
    """
    # 상관계수 행렬 cv (A, B)
    cv = np.array([
        [0.0, 0.0], [0.0, 0.0],              # F2, Cl2 (데이터 없음)
        [-71.9000, 0.2500],                  # SO2 (인덱스 3)
        [-26.5000, -21.3000], [-94.2000, -0.4200], [-22.3000, -1.7200],
        [-12.3000, 27.2000], [-58.6000, 12.7000], [-33.2000, 26.4000],
        [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], # H2, N2, O2, C2H4
        [-20.1000, 24.9000], [-23.3000, 49.7000], [-28.8000, 74.7000],
        [16.7000, 45.9000], [7.8000, 68.7000], [18.5000, 70.6000],
        [-25.0000, 56.5000], [10.3000, 47.7000], [-35.1000, 139.0000],
        [24.2000, 38.6000], [-50.2000, 36.5000], [-24.8000, 26.9000],
        [-22.6000, 32.1000]
    ])

    # 계수 조정 (B값에 1e-3 곱함)
    wv = cv.copy()
    wv[:, 1] = cv[:, 1] * 1e-3

    # 화합물 ID 확인
    ind = compID(cname)
    if ind is None:
        print("Error: 화합물을 찾을 수 없습니다.")
        return None

    # T를 numpy 배열로 변환 (반복문 및 조건 처리를 위해)
    T = np.atleast_1d(T)
    gf_list = []

    # 깁스 자유 에너지 계산 로직
    # 원본 MATLAB 코드의 if ind ~= 3 로직 구현
    if ind != 3:
        # 일반 화합물: gf = A + B*T
        row = wv[ind - 1]
        gf = row[0] + row[1] * T
        return gf if gf.size > 1 else gf[0]
    else:
        # 이산화황(SO2) 특수 케이스 처리
        for t_val in T:
            if 298 <= t_val <= 717:
                # 298K~717K 범위: 원본 코드에서 ind=1의 계수를 사용하는 로직 유지
                # (참고: MATLAB 원본에서 ind=1은 cv의 첫 번째 행임)
                row_so2 = wv[0] 
                gfv = row_so2[0] + row_so2[1] * t_val
            else:
                # 717K~1500K 범위: 특정 계수 [-86.8, 17.7*1e-3] 사용
                gfv = -86.8 + (17.7 * 1e-3) * t_val
            gf_list.append(gfv)
        
        gf = np.array(gf_list)
        return gf if gf.size > 1 else gf[0]

# 사용 예시:
# print(gfree(500, 'CO2'))
# print(gfree([400, 800], 'SO2'))