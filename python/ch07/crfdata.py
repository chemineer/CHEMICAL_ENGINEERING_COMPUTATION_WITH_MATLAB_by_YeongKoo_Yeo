import numpy as np

def crfdata():
    """MATLAB의 crfdata.m 스크립트를 파이썬 함수로 구현"""
    data = {}
    data['t'] = 0.00254       # membrane thickness (cm)
    data['Pm'] = np.array([50, 5]) * 1e-10 # permeability (cm^3*cm/(s*cm^2*cmHg))
    data['alpa'] = data['Pm'][0] / data['Pm'][1] # ratio of permeabilities
    data['ph'] = 80           # feed side pressure (cmHg)
    data['pl'] = 20           # permeate side pressure (cmHg)
    data['r'] = data['pl'] / data['ph'] # pressure ratio
    data['qf'] = 1e4          # feed rate (cm^3/s(STP))
    data['xf'] = 0.5          # Feed composition (mole fraction)
    data['theta'] = None      # stage-cut (None은 빈 배열 대용)
    data['xr'] = 0.25         # desired reject composition (mole fraction)
    
    return data