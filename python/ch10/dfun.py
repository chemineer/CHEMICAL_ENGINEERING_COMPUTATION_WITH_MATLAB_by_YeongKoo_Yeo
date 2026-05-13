import numpy as np

def dfun(x):
    """
    dfun: 목적 함수, 등식 제약 조건, 부등식 제약 조건의 그레이디언트를 계산하는 함수
    """
    # x가 리스트 형태일 경우를 대비해 numpy 배열로 변환
    x = np.array(x)
    
    # 목적 함수의 그레이디언트 (df)
    # MATLAB: df = [2 2*x(2)]
    df = np.array([2, 2 * x[1]])
    
    # 등식 제약 조건 h(x) = 0 의 그레이디언트 (dh)
    # MATLAB: dh = [2*x(1) 2*x(2)]
    dh = np.array([2 * x[0], 2 * x[1]])
    
    # 부등식 제약 조건 gi(x) >= 0 의 그레이디언트 (dg)
    # MATLAB: dg = [1 0; -1 0; 0 1; 0 -1]
    dg = np.array([
        [1, 0],
        [-1, 0],
        [0, 1],
        [0, -1]
    ])
    
    # 결과 결합 (dfv = [df' dh' dg'])
    # MATLAB의 df'와 dh'는 열 벡터이므로, 파이썬에서도 열 형태로 변환하여 결합합니다.
    # df.reshape(-1, 1)은 (2,) 배열을 (2, 1) 배열로 만듭니다.
    # dg'는 dg의 전치 행렬이므로 dg.T를 사용합니다.
    dfv = np.hstack([df.reshape(-1, 1), dh.reshape(-1, 1), dg.T])
    
    return dfv

# --- 테스트 코드 ---
if __name__ == "__main__":
    test_x = [1.0, 2.0]
    result = dfun(test_x)
    print("dfun result:")
    print(result)