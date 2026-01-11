import numpy as np
from engine.geometria import GeometriaTorre
from typing import Dict, Any

class ParametrosLinha:
    def __init__(self, geometria: GeometriaTorre, condutor_fase: Dict[str, Any], 
                 condutor_pr: Dict[str, Any], num_cond_fase: int = 1, 
                 espacamento_fase: float = 0.45, freq: float = 60.0, 
                 rho_solo: float = 100.0):
        self.geometria = geometria
        self.condutor_fase = condutor_fase
        self.condutor_pr = condutor_pr
        self.num_cond_fase = num_cond_fase
        self.espacamento_fase = espacamento_fase
        self.freq = freq
        self.rho_solo = rho_solo
        
        # Parâmetros de cálculo
        self.omega = 2 * np.pi * freq
        self.mu_0 = 4e-7 * np.pi  # H/m
        self.epsilon_0 = 8.8541878128e-12  # F/m
        
        # RMG e Raio Equivalente (convertendo mm para m)
        self.rmg_eq_fase = self.geometria.calcular_rmg_equivalente(
            condutor_fase['rmg_mm'] / 1000, 
            condutor_fase['diametro_mm'] / 2000,
            num_cond_fase, espacamento_fase
        )
        self.raio_eq_fase = self.geometria.calcular_raio_equivalente_bundle(
            condutor_fase['diametro_mm'] / 2000,
            num_cond_fase, espacamento_fase
        )
        
        self.rmg_pr = condutor_pr['rmg_mm'] / 1000
        self.raio_pr = condutor_pr['diametro_mm'] / 2000
        
        # Resistências (Ω/km)
        # Note: In a real simulation, we'd adjust for temperature as per planning
        self.r_fase = (condutor_fase['rac_25C'] / num_cond_fase) 
        self.r_pr = condutor_pr['rac_25C']

    def matriz_impedancia_primitiva(self) -> np.ndarray:
        """
        Calcula a matriz de impedância primitiva (Z_prim) usando o método de Dubanton (Ω/km)
        """
        coords = self.geometria.obter_coordenadas_completas()
        n = len(coords)
        
        # Profundidade de retorno de Dubanton
        d_e = 658.5 * np.sqrt(self.rho_solo / self.freq)
        
        z_prim = np.zeros((n, n), dtype=complex)
        
        # RMGs e Resistências para cada condutor na matriz
        # Ordem de coords: fases C1, [fases C2], para-raios
        n_fases_total = len(self.geometria.fases)
        if self.geometria.tipo_circuito == 'duplo':
            n_fases_total += len(self.geometria.fases_c2)
            
        rmgs = [self.rmg_eq_fase] * n_fases_total + [self.rmg_pr] * len(self.geometria.para_raios)
        resistencias = [self.r_fase] * n_fases_total + [self.r_pr] * len(self.geometria.para_raios)
        
        for i in range(n):
            xi, yi = coords[i]
            for j in range(n):
                xj, yj = coords[j]
                
                if i == j:
                    # Autoimpedância: Ri + j*ω*μ0/2π * ln(De / RMG_i)
                    # 1e3 factor for Ω/km
                    z_prim[i, i] = resistencias[i] + 1j * self.omega * (self.mu_0 / (2 * np.pi)) * np.log(d_e / rmgs[i]) * 1e3
                else:
                    # Impedância mútua: j*ω*μ0/2π * ln(Dij_img / Dij)
                    d_ij = np.sqrt((xi - xj)**2 + (yi - yj)**2)
                    # Imagem j está em (xj, -yj - 2*p) onde p é profundidade complexa, 
                    # mas Dubanton simplifica Dij_img ≈ sqrt((xi-xj)^2 + (yi+yj+2*p)^2)
                    # Para Dubanton real: De ≈ |yi + yj + 2*p|
                    d_ij_img = np.sqrt((xi - xj)**2 + (yi + yj + 2 * (d_e / np.sqrt(2)))**2) # Aproximação comum
                    # Uma forma mais rigorosa de Dubanton: Dij_img = sqrt((xi-xj)^2 + (yi+yj+2*p)^2)
                    # Onde p = sqrt(rho / (j * omega * mu0))
                    z_prim[i, j] = 1j * self.omega * (self.mu_0 / (2 * np.pi)) * np.log(d_ij_img / d_ij) * 1e3
                    
        return z_prim

    def matriz_potencial_maxwell(self) -> np.ndarray:
        """
        Calcula a matriz de coeficientes de potencial de Maxwell (km/F)
        """
        coords = self.geometria.obter_coordenadas_completas()
        n = len(coords)
        
        n_fases_total = len(self.geometria.fases)
        if self.geometria.tipo_circuito == 'duplo':
            n_fases_total += len(self.geometria.fases_c2)
            
        raios = [self.raio_eq_fase] * n_fases_total + [self.raio_pr] * len(self.geometria.para_raios)
        
        p = np.zeros((n, n))
        
        factor = 1 / (2 * np.pi * self.epsilon_0)
        
        for i in range(n):
            xi, yi = coords[i]
            for j in range(n):
                xj, yj = coords[j]
                
                if i == j:
                    # Pii = 1/2πε * ln(2h/r)
                    p[i, i] = factor * np.log(2 * yi / raios[i])
                else:
                    # Pij = 1/2πε * ln(Dij_img / Dij)
                    d_ij = np.sqrt((xi - xj)**2 + (yi - yj)**2)
                    d_ij_img = np.sqrt((xi - xj)**2 + (yi + yj)**2)
                    p[i, j] = factor * np.log(d_ij_img / d_ij)
                    
        return p * 1e-3 # Convert to km/F (so C will be F/km)

    def reducao_kron(self, matriz: np.ndarray) -> np.ndarray:
        """
        Elimina os para-raios da matriz primitava para obter a matriz reduzida (fases)
        """
        n_fases = len(self.geometria.fases)
        if self.geometria.tipo_circuito == 'duplo':
            n_fases += len(self.geometria.fases_c2)
            
        z_ff = matriz[:n_fases, :n_fases]
        z_fg = matriz[:n_fases, n_fases:]
        z_gf = matriz[n_fases:, :n_fases]
        z_gg = matriz[n_fases:, n_fases:]
        
        if len(z_gg) == 0:
            return z_ff
            
        z_abc = z_ff - z_fg @ np.linalg.inv(z_gg) @ z_gf
        return z_abc

    def calcular_admitancia(self) -> np.ndarray:
        """
        Retorna a matriz de admitância reduzida Y_abc (S/km)
        Y = G + jB ≈ jωC
        """
        p_prim = self.matriz_potencial_maxwell()
        c_prim = np.linalg.inv(p_prim)
        c_abc = self.reducao_kron(c_prim)
        
        # Admitância capacitiva: j * omega * C
        y_abc = 1j * self.omega * c_abc
        return y_abc
