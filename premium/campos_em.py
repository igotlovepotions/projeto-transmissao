import numpy as np
from typing import Dict, List, Tuple

class CamposEM:
    def __init__(self, coords: List[Tuple[float, float]], tensoes_fn: List[complex], correntes: List[complex]):
        """
        coords: Lista de (x, y) dos condutores
        tensoes_fn: Lista de tensões fase-neutro (kV)
        correntes: Lista de correntes (A)
        """
        self.coords = coords
        self.tensoes = tensoes_fn
        self.correntes = correntes
        self.epsilon_0 = 8.854e-12
        self.mu_0 = 4e-7 * np.pi

    def campo_eletrico_solo(self, x_ponto: float) -> float:
        """
        Calcula o campo elétrico no solo (y=0) usando o método das imagens.
        Retorna valor eficaz (kV/m).
        """
        e_total = 0j
        
        for i, (xc, yc) in enumerate(self.coords):
            # Carga linear aproximada q = 2πε * V / ln(2h/r)
            # Para simplificar, vamos usar uma formulação direta de potencial
            # Er = V / (r * ln(H/r)) - Simplificado
            
            # Formulação de Maxwell: Sum(Qi * ln(Dij_img/Dij)) = 2πε * Vi
            # Aqui simplificamos assumindo que conhecemos as tensões
            
            r1 = np.sqrt((x_ponto - xc)**2 + yc**2)
            r2 = np.sqrt((x_ponto - xc)**2 + yc**2) # Imagem está na mesma distância de (x, 0)
            
            # Contribuição vertical (Ey) é a que resta no solo
            # Ey = (q / 2πε) * [ yc / r1^2 + yc / r2^2 ]
            # Usando V/ln(...) como proxy para q/2πε
            q_factor = self.tensoes[i] / 5.0 # Fator de escala aproximado para exemplo
            
            e_cont = q_factor * (yc / (r1**2) + yc / (r2**2))
            e_total += e_cont
            
        return abs(e_total)

    def campo_magnetico_solo(self, x_ponto: float) -> float:
        """
        Lei de Biot-Savart para condutores infinitos.
        Retorna densidade de fluxo B (uT).
        """
        b_total_x = 0j
        b_total_y = 0j
        
        for i, (xc, yc) in enumerate(self.coords):
            r2 = (x_ponto - xc)**2 + yc**2
            
            # B = mu0 * I / (2pi * r)
            # Bx = -B * sin(theta) = -B * (yc/r)
            # By = B * cos(theta) = B * (x-xc)/r
            
            factor = (self.mu_0 * self.correntes[i]) / (2 * np.pi * r2)
            
            b_total_x += -factor * yc
            b_total_y += factor * (x_ponto - xc)
            
        b_res = np.sqrt(abs(b_total_x)**2 + abs(b_total_y)**2)
        return b_res * 1e6 # Tesla -> microTesla
