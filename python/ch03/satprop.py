from satH2Oprop import satH2Oprop

def satprop():
    """
    포화 증기/액체 상태의 H2O 물성을 출력하는 메인 스크립트
    """
    
    # 사용자로부터 입력 받기
    # MATLAB의 input()은 변수를 바로 인식하지만, 파이썬은 문자열 입력 시 따옴표가 필요합니다.
    dtype = input("Data type ('P' for Pressure in kPa, 'T' for Temperature in deg.C) = ")
    dvalue = float(input("Data value = "))
    phH2O = input("Phase ('V' for Vapor, 'L' for Liquid) = ")
    
    # satH2Oprop 함수는 별도로 정의되어 있어야 합니다.
    # 해당 함수가 딕셔너리 형태의 결과값을 반환한다고 가정합니다.
    w = satH2Oprop(dtype, phH2O, dvalue)
    
    # 결과 출력
    print(f"\nSaturation pressure = {w['P']} kPa")
    print(f"Saturation temperature = {w['T']} deg.C")
    print(f"Specific enthalpy = {w['H']} kJ/kg")
    print(f"Specific entropy = {w['S']} kJ/(kg-K)")
    print(f"Specific internal energy = {w['U']} kJ/kg")
    print(f"Specific volume = {w['V']} cm^3/g")

# 사용 방법:
# satprop()