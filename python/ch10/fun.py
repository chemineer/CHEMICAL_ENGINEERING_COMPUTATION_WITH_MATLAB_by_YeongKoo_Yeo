import numpy as np

def fun(x):
    """
    fun: 목적 함수 및 제약 조건 계산 함수
    입력 x는 [x[0], x[1]] 형태의 배열 또는 리스트여야 합니다.
    """
    x = np.array(x)
    
    # 각 요소 계산
    # 1. 목적 함수
    obj = 2 * x[0] + x[1]**2
    
    # 2. 등식 제약 조건: h(x) = 0
    h = x[0]**2 + x[1]**2 - 8
    
    # 3. 부등식 제약 조건: g1(x) >= 0, g2(x) >= 0, g3(x) >= 0, g4(x) >= 0
    g1 = x[0]
    g2 = -x[0] + 4
    g3 = x[1] - 1
    g4 = -x[1] + 5
    
    # 모든 결과를 리스트(또는 numpy 배열)로 결합하여 반환
    fv = np.array([obj, h, g1, g2, g3, g4])
    
    return fv

# --- 테스트 코드 ---
if __name__ == "__main__":
    test_x = [1.0, 2.0]
    result = fun(test_x)
    print("fun result:")
    print(result)