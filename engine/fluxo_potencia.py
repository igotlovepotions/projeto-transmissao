import numpy as np
from typing import Dict, Any, Tuple

class FluxoPotencia:
    def __init__(self, abcd_mat: np.ndarray, v_nominal_kv: float):
        self.abcd = abcd_mat
        self.v_nom = v_nominal_kv / np.sqrt(3) # V_ln (kV)

    def calcular_terminal_receptor(self, vs_pu: complex, p_rec_mw: float, fp: float, tipo_fp: str = 'atraso') -> Dict[str, Any]:
        """
        Dado Vs (pu), calcula Vr e Ir no terminal receptor para uma carga Pr.
        Considera sistema trifásico, cálculos em pu ou valores reais por fase.
        """
        a, b, c, d = self.abcd[0,0], self.abcd[0,1], self.abcd[1,0], self.abcd[1,1]
        
        # Tensão de envio em kV (fase-neutro)
        vs = vs_pu * self.v_nom 
        
        # Potência recebida por fase (MW)
        pr_fase = p_rec_mw / 3
        
        # Iteração para encontrar Vr
        # Vs = A*Vr + B*Ir
        # Ir = (Pr - jQr)* / Vr*
        
        phi = np.arccos(fp)
        if tipo_fp == 'atraso':
            # Qr positivo (indutivo consome reativos)
            qr_fase = pr_fase * np.tan(phi)
        else:
            # Qr negativo (capacitivo fornece reativos)
            qr_fase = -pr_fase * np.tan(phi)
            
        vr = vs / abs(a) # Chute inicial
        
        for _ in range(20):
            sr = pr_fase + 1j * qr_fase
            ir = np.conj(sr / vr)
            vr_novo = (vs - b * ir) / a
            
            if abs(vr_novo - vr) < 1e-6:
                vr = vr_novo
                break
            vr = vr_novo
            
        # Corrente de envio
        is_src = c * vr + d * ir
        ss_src = vs * np.conj(is_src)
        
        return {
            'Vr_kv': vr,
            'Vr_pu': abs(vr) / self.v_nom,
            'Ir_a': ir * 1e3, # A
            'Is_a': is_src * 1e3,
            'Pr_mw': pr_fase * 3,
            'Ps_mw': ss_src.real * 3,
            'Perdas_mw': (ss_src.real - pr_fase) * 3,
            'Eficiencia': (pr_fase / ss_src.real) * 100,
            'Regulacao': ((abs(vs/a) - abs(vr)) / abs(vr)) * 100
        }

    def calcular_sil(self, zc: complex) -> float:
        """Surge Impedance Loading (MW)"""
        # SIL = Vnom^2 / |Zc|
        # Vnom em kV (linha-linha)
        v_ll = self.v_nom * np.sqrt(3)
        sil = (v_ll**2) / abs(zc)
        return sil

    def curva_potencia_angulo(self, vs_kv_ln: float, vr_kv_ln: float) -> Tuple[np.ndarray, np.ndarray]:
        """Dados Vs e Vr, retorna curva P x delta"""
        delta = np.linspace(0, 180, 181)
        a_abs = abs(self.abcd[0,0])
        a_ang = np.angle(self.abcd[0,0])
        b_abs = abs(self.abcd[0,1])
        b_ang = np.angle(self.abcd[0,1])
        
        # Pr = (|Vs||Vr|/|B|)cos(beta-delta) - (|A||Vr|^2/|B|)cos(beta-alpha)
        # Onde beta = ang(B), alpha = ang(A)
        
        delta_rad = np.radians(delta)
        
        p = ( (vs_kv_ln * vr_kv_ln / b_abs) * np.cos(b_ang - delta_rad) - 
              (a_abs * vr_kv_ln**2 / b_abs) * np.cos(b_ang - a_ang) )
        
        return delta, p * 3 # Total trifásico (MW)
