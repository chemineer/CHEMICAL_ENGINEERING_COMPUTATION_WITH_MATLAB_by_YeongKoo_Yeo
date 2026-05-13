def compID(cname):
    """
    화합물 이름을 입력받아 해당 ID 번호를 반환하는 함수
    """
    # 입력값을 대문자로 변환
    cname = cname.upper()
    
    # 화합물과 ID 매핑 딕셔너리
    comp_map = {
        'FLUORINE': 1, 'F2': 1,
        'CHLORINE': 2, 'CL2': 2,
        'SULFUR DIOXIDE': 3, 'SO2': 3,
        'CARBON MONOXIDE': 4, 'CO': 4,
        'CARBON DIOXIDE': 5, 'CO2': 5,
        'HYDROGEN CHLORIDE': 6, 'HCL': 6,
        'AMMONIA': 7, 'NH3': 7,
        'WATER': 8, 'H2O': 8,
        'HYDROGEN PEROXIDE': 9, 'H2O2': 9,
        'HYDROGEN': 10, 'H2': 10,
        'NITROGEN': 11, 'N2': 11,
        'OXYGEN': 12, 'O2': 12,
        'ETHYLENE': 13, 'C2H4': 13,
        'METHANE': 14, 'CH4': 14,
        'ETHANE': 15, 'C2H6': 15,
        'PROPANE': 16, 'C3H8': 16,
        'BENZENE': 17, 'C6H6': 17,
        'TOLUENE': 18, 'C7H8': 18,
        'ANILINE': 19, 'C6H7N': 19,
        'PHENOL': 20, 'C6H6O': 20,
        'CYCLOPROPANE': 21, 'C3H6': 21,
        'CYCLOHEXANE': 22, 'C6H12': 22,
        '1,3 BUTADIENE': 23, 'C4H6': 23,
        'METHANOL': 24, 'CH4O': 24,
        'CHLOROFORM': 25, 'CHCL3': 25,
        'CARBON TETRACHLORIDE': 26, 'CCL4': 26
    }
    
    # 해당 화합물이 딕셔너리에 있으면 ID 반환, 없으면 None 반환
    return comp_map.get(cname)

# 사용 예시:
# print(compID('CO2'))  # 출력: 5