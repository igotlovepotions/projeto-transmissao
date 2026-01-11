import numpy as np
from typing import Dict, Any

class EfeitoCorona:
    def __init__(self, raio_cm: float, dmd_cm: float):
        self.r = raio_cm
        self.d = dmd_cm

    def calcular_tensao_critica_peek(self, altitude_m: int, condicao_superficie: str = 'sujo') -> float:
        """
        Fórmula de Peek para tensão crítica disruptiva (kV fase-neutro).
        """
        # Fator de rugosidade m0
        m0_map = {
            'liso': 1.0,
            'oxidado': 0.98,
            'sujo': 0.82 # Valor conservador para cabos encordoados
        }
        m0 = m0_map.get(condicao_superficie, 0.82)
        
        # Densidade relativa do ar (delta)
        delta = np.exp(-altitude_m / 8150)
        
        # Gradiente crítico (kV/cm)
        # E0 = 21.1 * m0 * delta * (1 + 0.3 / sqrt(delta * r))
        e0 = 21.1 * m0 * delta * (1 + 0.3 / np.sqrt(delta * self.r))
        
        # Tensão crítica (kV rms fn)
        v_crit = e0 * self.r * np.log(self.d / self.r)
        
        return v_crit

    def perdas_fowler(self, v_op_kv: float, v_crit_kv: float, freq: float = 60.0) -> float:
        """Estimativa de perdas por corona (kW/km)"""
        if v_op_kv <= v_crit_kv:
            return 0.0
        
        # Fórmula simplificada: P = k * f * (V - Vc)^2
        return 0.05 * freq * (v_op_kv - v_crit_kv)**2
