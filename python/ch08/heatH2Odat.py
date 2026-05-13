import math

def heatH2Odat():
    """
    heatH2Odat.m: Physical properties and operating conditions for H2O heater
    Returns a dictionary containing all parameters.
    """
    data = {}

    # Physical properties (hot(shell): demineralized water, cold(tube): raw water)
    data['Cph'] = 1.0
    data['Cpc'] = 1.01
    data['rhoh'] = 62.4
    data['rhoc'] = 62.4
    data['kh'] = 0.36
    data['kc'] = 0.363
    data['muh'] = 0.81
    data['muc'] = 0.92
    data['muwc'] = data['muc']
    data['muwh'] = data['muh']

    # Operating conditions
    data['mh'] = 1.5e5
    data['ui'] = 5
    data['T1'] = 95
    data['T2'] = 85
    data['t1'] = 75
    data['t2'] = 80

    # Shell and tube dimensions
    data['L'] = 10
    data['Do'] = 0.0625
    data['Di'] = 0.05167
    data['Ds'] = 1.77083
    data['Pt'] = 1/12
    data['cl'] = 0.02083
    data['Rfi'] = 0.001
    data['Rfo'] = 0
    data['kw'] = 30

    # Parameters and basic properties
    data['gc'] = 32.2
    
    # Heat load (Btu/h)
    data['Qr'] = data['mh'] * data['Cph'] * (data['T1'] - data['T2'])
    
    # Cold stream rate (lb/h)
    data['mc'] = data['Qr'] / (data['Cpc'] * (data['t2'] - data['t1']))
    
    # Log-mean temperature difference (dTlm)
    dT1 = data['T1'] - data['t2']
    dT2 = data['T2'] - data['t1']
    data['dTlm'] = (dT1 - dT2) / math.log(dT1 / dT2)
    
    # Correction factor (F12) calculation
    R = (data['T1'] - data['T2']) / (data['t2'] - data['t1'])
    S = (data['t2'] - data['t1']) / (data['T1'] - data['t1'])
    data['R'] = R
    data['S'] = S
    
    sqrt_R2_1 = math.sqrt(R**2 + 1)
    F12den = (R - 1) * math.log((2 - S * (R + 1 - sqrt_R2_1)) / (2 - S * (R + 1 + sqrt_R2_1)))
    data['F12'] = sqrt_R2_1 * math.log((1 - S) / (1 - R * S)) / F12den
    
    # Tube-side total cross-sectional area (ft^2/pass)
    data['Aci'] = data['mc'] / (data['rhoc'] * data['ui'] * 3600)
    
    # Number of tubes per pass
    data['Nt'] = math.ceil(4 * data['Aci'] / (math.pi * data['Di']**2))
    
    # Heat transfer area per tube
    data['At'] = math.pi * data['Di'] * data['L']
    
    # Tube thickness and diameters
    data['dx'] = (data['Do'] - data['Di']) / 2
    data['Dlm'] = (data['Do'] - data['Di']) / math.log(data['Do'] / data['Di'])
    
    # Hydraulic effective diameter (De)
    data['De'] = (4 / (math.pi * data['Do'])) * (data['Pt']**2 - math.pi * data['Do']**2 / 4)
    
    # Initial assumed Overall Heat Transfer Coefficient
    data['Ui_init'] = 400 
    
    return data

# 사용 예시:
# params = heatH2Odat()
# print(params['Qr'])