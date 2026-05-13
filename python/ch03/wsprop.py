# wsprop.py: H2O 속성(H, S, U, V) 계산 스크립트

def wsprop():
    # 사용자로부터 입력 받기
    try:
        P = float(input('Pressure (kPa) = '))
        T = float(input('Temperature (deg.C) = '))
        
        # H2Oprop 함수 호출 (별도의 H2Oprop.py가 정의되어 있어야 함)
        # 딕셔너리 형태로 속성값들을 반환받는다고 가정합니다.
        w = H2Oprop(P, T)
        
        # 결과 출력
        print(f'Specific enthalpy = {w["H"]} kJ/kg')
        print(f'Specific entropy = {w["S"]} kJ/(kg-K)')
        print(f'Specific internal energy = {w["U"]} kJ/kg')
        print(f'Specific volume = {w["V"]} cm^3/g')
        
    except ValueError:
        print("유효한 숫자를 입력해주세요.")

# 참고: H2Oprop 함수는 별도의 파일로 존재해야 합니다.
# 예시 구조:
# def H2Oprop(P, T):
#     # 여기에 계산 로직 구현
#     return {'H': ..., 'S': ..., 'U': ..., 'V': ...}

if __name__ == "__main__":
    wsprop()