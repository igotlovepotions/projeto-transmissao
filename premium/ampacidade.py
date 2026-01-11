import numpy as np
from typing import Dict, Any

class AmpacidadeLinha:
    def __init__(self, diametro_mm: float, rac_temp_max: float, temp_max_c: float = 75.0):
        self.d_mm = diametro_mm
        self.r_max = rac_temp_max # Resistência na temp_max (Ω/km)
        self.temp_max = temp_max_c
        
    def calcular_limite_termico(self, temp_amb_c: float, vento_m_s: float, radiacao_w_m2: float = 1000.0) -> Dict[str, Any]:
        """
        Balanço térmico IEEE 738:
        I^2 * R = Qc + Qr - Qs
        Qc: Convecção, Qr: Radiação, Qs: Ganho Solar
        """
        d_m = self.d_mm / 1000
        delta_t = self.temp_max - temp_amb_c
        
        if delta_t <= 0:
            return {'I_max_A': 0.0, 'status': 'Temp Amb > Temp Max'}

        # 1. Perdas por Convecção (Simplificado)
        # Qc = [ 1.01 + 1.35 * (Re^0.52) ] * k * delta_t
        # Usando fórmula prática comum:
        h_conv = 10.0 + 6.0 * np.sqrt(vento_m_s) # Coeficiente aprox (W/m2.K)
        qc = h_conv * (np.pi * d_m) * delta_t
        
        # 2. Perdas por Radiação
        sigma = 5.67e-8
        emissividade = 0.5
        t_abs_max = self.temp_max + 273.15
        t_abs_amb = temp_amb_c + 273.15
        qr = sigma * emissividade * (np.pi * d_m) * (t_abs_max**4 - t_abs_amb**4)
        
        # 3. Ganho Solar
        absortividade = 0.5
        qs = absortividade * d_m * radiacao_w_m2
        
        # Balanço de potência por metro: I^2 * (R/1000) = qc + qr - qs
        r_metro = self.r_max / 1000
        
        perdas_termicas = qc + qr - qs
        if perdas_termicas < 0:
            return {'I_max_A': 0.0, 'status': 'Ganho solar excessivo'}
            
        i_max = np.sqrt(perdas_termicas / r_metro)
        
        return {
            'I_max_A': i_max,
            'Qc_w_m': qc,
            'Qr_w_m': qr,
            'Qs_w_m': qs,
            'R_max_ohm_m': r_metro
        }
