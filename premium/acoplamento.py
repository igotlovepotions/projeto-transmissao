import numpy as np
from typing import Dict, Any

class AcoplamentoMutuo:
    def __init__(self, z_prim_6x6: np.ndarray):
        """
        z_prim_6x6: Matriz de impedância primitiva para dois circuitos trifásicos (6x6) 
        após eliminação dos para-raios por Kron.
        """
        self.z_6x6 = z_prim_6x6

    def extrair_blocos(self) -> Dict[str, np.ndarray]:
        """
        Extrai as submatrizes de autoimpedância de cada circuito e a mútua entre eles.
        """
        z11 = self.z_6x6[0:3, 0:3]
        z22 = self.z_6x6[3:6, 3:6]
        z12 = self.z_6x6[0:3, 3:6]
        
        return {
            'Z11': z11,
            'Z22': z22,
            'Z12': z12
        }

    def calcular_sequencia_mutua(self) -> complex:
        """
        Calcula a impedância de sequência zero mútua (Z0m/Z012_mutua).
        Ponto crítico para estudos de proteção e faltas em circuitos duplos.
        """
        a = np.exp(1j * 2 * np.pi / 3)
        t_inv = (1/3) * np.array([
            [1, 1, 1],
            [1, a, a**2],
            [1, a**2, a]
        ])
        t = np.array([
            [1, 1, 1],
            [1, a**2, a],
            [1, a, a**2]
        ])
        
        z12 = self.z_6x6[0:3, 3:6]
        z012_mutua = t_inv @ z12 @ t
        
        return z012_mutua[0, 0] # Retorna sequência zero mútua
