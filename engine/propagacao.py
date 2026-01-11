import numpy as np
from typing import Dict, Union

class PropagacaoLinha:
    def __init__(self, z_seq: Dict[str, complex], y_seq: Dict[str, complex], freq: float = 60.0):
        self.z_seq = z_seq
        self.y_seq = y_seq
        self.freq = freq
        
    def calcular_impedancia_caracteristica(self) -> Dict[str, complex]:
        """Zc = sqrt(Z / Y)"""
        zc_seq = {}
        for seq in ['0', '1', '2']:
            if self.y_seq[seq] == 0:
                zc_seq[seq] = float('inf')
            else:
                zc_seq[seq] = np.sqrt(self.z_seq[seq] / self.y_seq[seq])
        return zc_seq

    def calcular_constante_propagacao(self) -> Dict[str, Dict[str, Union[complex, float]]]:
        """gamma = sqrt(Z * Y) = alpha + j*beta"""
        gamma_seq = {}
        for seq in ['0', '1', '2']:
            gamma = np.sqrt(self.z_seq[seq] * self.y_seq[seq])
            gamma_seq[seq] = {
                'gamma': gamma,
                'alpha': gamma.real,  # Atenuação (Np/km)
                'beta': gamma.imag    # Constante de fase (rad/km)
            }
        return gamma_seq

    def calcular_parametros_onda(self, beta: float) -> Dict[str, float]:
        """v = omega / beta, lambda = 2pi / beta"""
        omega = 2 * np.pi * self.freq
        if beta == 0:
            return {'velocidade_km_s': 0, 'velocidade_pct_luz': 0, 'comprimento_onda_km': 0}
            
        v = omega / beta
        v_pct_luz = (v / 299792.458) * 100 # v / c * 100
        lambda_onda = 2 * np.pi / beta
        
        return {
            'velocidade_km_s': v,
            'velocidade_pct_luz': v_pct_luz,
            'comprimento_onda_km': lambda_onda
        }

    def matriz_abcd(self, seq: str, comprimento_km: float, modelo: str = 'auto') -> np.ndarray:
        """
        Calcula a matriz ABCD para uma determinada sequência e comprimento.
        modelos: 'curta' (< 80km), 'media' (80-250km), 'longa' (> 250km)
        """
        l = comprimento_km
        z_unit = self.z_seq[seq]
        y_unit = self.y_seq[seq]
        
        if modelo == 'auto':
            if l < 80: modelo = 'curta'
            elif l < 250: modelo = 'media'
            else: modelo = 'longa'

        if modelo == 'curta':
            a = d = 1
            b = z_unit * l
            c = 0
        elif modelo == 'media':
            # Modelo Pi Nominal
            z_total = z_unit * l
            y_total = y_unit * l
            a = d = 1 + (z_total * y_total) / 2
            b = z_total
            c = y_total * (1 + (z_total * y_total) / 4)
        else:
            # Modelo de Linha Longa (Equações de Propagação)
            gamma = np.sqrt(z_unit * y_unit)
            zc = np.sqrt(z_unit / y_unit)
            
            a = d = np.cosh(gamma * l)
            b = zc * np.sinh(gamma * l)
            c = np.sinh(gamma * l) / zc
            
        return np.array([[a, b], [c, d]])
