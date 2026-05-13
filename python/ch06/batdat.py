def batdat():
    """
    배치 반응기(batch reactor) 시뮬레이션 데이터 반환 함수
    
    Returns:
    data: 반응기 파라미터가 담긴 딕셔너리
    """
    data = {
        "A1": 1.2,
        "A2": 180.0,
        "E1": 2.1e4,
        "E2": 4.3e4,
        "R": 8.314,
        "dH1": 4.09e4,
        "dH2": 8.24e4,
        "rho": 1000.0,
        "Cp": 1.0,
        "Tc": 20.0,
        "Ts": 110.0,
        "Uj": 1.2,
        "Uc": 3.0,
        "AcV": 18.6,
        "AjV": 31.5
    }
    return data

# 사용 예시:
# params = batdat()
# print(params["A1"])