import numpy as np

def grgfun(x):
    """
    grgfun: 목적 함수(f)와 제약 조건(g)을 계산하는 함수
    x1~x4: 결정 변수, x5~x7: 슬랙(Slack) 변수
    """
    # x가 리스트나 튜플일 경우를 대비해 numpy 배열로 변환
    x = np.array(x)
    
    # MATLAB 인덱스 1~7을 파이썬 인덱스 0~6으로 매핑
    x1, x2, x3, x4 = x[0], x[1], x[2], x[3]
    x5, x6, x7 = x[4], x[5], x[6]
    
    # 목적 함수 f 계산
    f = (x1**2 + x2**2 + 2.3 * x3**2 - 1.2 * x4**2 - 
         4*x1 - 6*x2 - 20*x3 + 6*x4 + 100)
    
    # 제약 조건 g 계산 (결과를 담을 리스트 또는 배열 생성)
    g = np.zeros(3)
    
    # g(1)
    g[0] = (x1**2 + x2**2 + x3**2 + x4**2 + 
            x1 - x2 + x3 - x4 + x5 - 7)
    
    # g(2)
    g[1] = (x1**2 + 2 * x2**2 + x3**2 + 2 * x4**2 - 
            x1 - x4 + x6 - 11)
    
    # g(3)
    g[2] = (2 * x1**2 + x2**2 + x3**2 + 
            2 * x1 - x2 - x4 + x7 - 6)
    
    return f, g

# --- 사용 예시 ---
if __name__ == "__main__":
    # 7개의 변수를 가진 입력 벡터 예시
    test_x = [1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
    f_val, g_vals = grgfun(test_x)
    
    print(f"Objective function (f): {f_val}")
    print(f"Constraints (g): {g_vals}")