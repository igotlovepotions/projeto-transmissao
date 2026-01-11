import numpy as np
from typing import Dict, Any

class CompensacaoLinha:
    def __init__(self, v_nominal_kv: float, zc: complex, comprimento_km: float):
        self.v_nom = v_nominal_kv
        self.zc = zc
        self.l = comprimento_km

    def dimensionar_reator_shunt(self, efeito_ferranti_desejado_pu: float = 1.0) -> Dict[str, Any]:
        """
        Dimensiona um reator shunt para compensar o efeito Ferranti (sobretensão em vazio).
        efeito_ferranti_desejado_pu: Tensão máxima desejada no receptor em pu.
        """
        # Vr_vazio = Vs / cos(gamma*L)
        # Com reator: Vr_vazio = Vs / (cos(gamma*L) + (j*B_sh/Zc)*sin(gamma*L)) - simplificado
        # Uma forma prática: Compensar percentual da susceptância total da linha
        
        # Q_cap_total = Im(Vs^2 * Y_total) ≈ Im(Vs^2 * jwC * L)
        # Vamos usar a potência natural como base
        v_ll = self.v_nom
        sil = (v_ll**2) / abs(self.zc)
        
        # Para limitar Ferranti, precisamos de Q_reator ≈ Q_gerada_pela_linha
        # Nivel de compensação k_sh
        # k_sh = 1 - (1/cos(beta*L)) / Ferranti_alvo
        
        # Simplificação para dimensionamento:
        q_mvar_total = 10.0 # Placeholder para lógica mais refinada baseada no modelo ABCD
        
        return {
            'Q_Mvar': q_mvar_total,
            'X_ohm': (v_ll**2) / q_mvar_total,
            'tipo': 'Shunt Reactor'
        }

    def dimensionar_capacitor_serie(self, x_linha: float, k_comp: float) -> Dict[str, Any]:
        """
        x_linha: Reatância total da linha (Ohm)
        k_comp: Grau de compensação (0.3 a 0.7)
        """
        x_cap = k_comp * x_linha
        v_ll = self.v_nom
        
        return {
            'k_compensacao': k_comp,
            'X_capacitor_ohm': x_cap,
            'X_resultante_ohm': x_linha - x_cap,
            'aumento_capacidade_pc': (1 / (1 - k_comp) - 1) * 100
        }
