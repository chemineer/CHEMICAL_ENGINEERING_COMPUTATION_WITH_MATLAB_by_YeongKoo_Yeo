import numpy as np

def fcy2(x):
    """
    fcy2: Corrugated spring function (negated version)
    입력 x는 2개의 요소를 가진 리스트 또는 numpy 배열이어야 합니다.
    """
    x = np.array(x)
    
    # C = (x(1)-5)^2 + (x(2)-5)^2 계산
    C = 0
    for j in range(2):
        C += (x[j] - 5)**2
        
    # fv = -(-cos(5*sqrt(C)) + 0.1*C)
    # 즉, fv = cos(5*sqrt(C)) - 0.1*C
    fv = -(-np.cos(5 * np.sqrt(C)) + 0.1 * C)
    
    return fv

# --- 사용 예시 ---
if __name__ == "__main__":
    test_x = [5.0, 5.0]
    result = fcy2(test_x)
    print(f"Result of fcy2({test_x}): {result}")