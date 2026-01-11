import numpy as np
import plotly.graph_objects as go
from typing import Dict, Tuple, List, Optional

class GeometriaTorre:
    def __init__(self, tipo_circuito: str = 'simples'):
        """
        tipo_circuito: 'simples' ou 'duplo'
        """
        self.tipo_circuito = tipo_circuito
        self.fases: Dict[str, Tuple[float, float]] = {}      # {'A': (x, y), 'B': (x, y), 'C': (x, y)}
        self.para_raios: Dict[str, Tuple[float, float]] = {} # {'PR1': (x, y), 'PR2': (x, y)}
        self.fases_c2: Dict[str, Tuple[float, float]] = {}   # Para circuito duplo: {'A2': (x, y), ...}
        
    def adicionar_fase(self, nome: str, x: float, y: float, circuito: int = 1):
        if circuito == 1:
            self.fases[nome] = (x, y)
        else:
            self.fases_c2[nome] = (x, y)
            
    def adicionar_para_raio(self, nome: str, x: float, y: float):
        self.para_raios[nome] = (x, y)

    def calcular_rmg_equivalente(self, rmg_condutor: float, raio_externo_condutor: float, num_condutores: int, espacamento: float) -> float:
        """
        Calcula o RMG equivalente para um feixe (bundle) de condutores.
        rmg_condutor: RMG do condutor individual (m)
        num_condutores: número de condutores no feixe (1 a 6)
        espacamento: distância entre condutores adjacentes (m)
        """
        if num_condutores == 1:
            return rmg_condutor
        
        # Raio do círculo do bundle (distância do centro ao centro de cada condutor)
        if num_condutores == 2:
            r_bundle = espacamento / 2
        elif num_condutores == 3:
            r_bundle = espacamento / (np.sqrt(3))
        elif num_condutores == 4:
            r_bundle = espacamento / (np.sqrt(2))
        elif num_condutores == 6:
            r_bundle = espacamento # Para hexágono regular
        else:
            # Aproximação para outros casos (polígono regular)
            r_bundle = espacamento / (2 * np.sin(np.pi / num_condutores))

        # RMG_eq = n√(n * RMG * r^(n-1)) - Versão Stevenson/Fuchs
        # Ou RMG_eq = n√(RMG * Produto_distancias)
        # Simplificação comum: RMG_eq = (n * rmg * r^(n-1))^(1/n)
        return (num_condutores * rmg_condutor * (r_bundle**(num_condutores-1)))**(1/num_condutores)

    def calcular_raio_equivalente_bundle(self, raio_condutor: float, num_condutores: int, espacamento: float) -> float:
        """
        Calcula o raio equivalente para capacitância (raio externo do bundle)
        """
        if num_condutores == 1:
            return raio_condutor
        
        if num_condutores == 2:
            r_bundle = espacamento / 2
        elif num_condutores == 3:
            r_bundle = espacamento / (np.sqrt(3))
        elif num_condutores == 4:
            r_bundle = espacamento / (np.sqrt(2))
        else:
            r_bundle = espacamento / (2 * np.sin(np.pi / num_condutores))
            
        return (num_condutores * raio_condutor * (r_bundle**(num_condutores-1)))**(1/num_condutores)

    def obter_coordenadas_completas(self) -> List[Tuple[float, float]]:
        """Retorna lista de (x, y) de todas as fases e para-raios na ordem: fases_c1, [fases_c2], para_raios"""
        coords = list(self.fases.values())
        if self.tipo_circuito == 'duplo':
            coords.extend(list(self.fases_c2.values()))
        coords.extend(list(self.para_raios.values()))
        return coords

    def plotar_torre_2d(self):
        fig = go.Figure()

        # Plot Fases Circuito 1
        x_f1, y_f1 = zip(*self.fases.values()) if self.fases else ([], [])
        fig.add_trace(go.Scatter(x=x_f1, y=y_f1, mode='markers+text', 
                                 name='Fases C1', text=list(self.fases.keys()),
                                 marker=dict(size=12, color='red')))

        # Plot Fases Circuito 2
        if self.tipo_circuito == 'duplo' and self.fases_c2:
            x_f2, y_f2 = zip(*self.fases_c2.values())
            fig.add_trace(go.Scatter(x=x_f2, y=y_f2, mode='markers+text', 
                                     name='Fases C2', text=list(self.fases_c2.keys()),
                                     marker=dict(size=12, color='orange')))

        # Plot Para-raios
        if self.para_raios:
            x_pr, y_pr = zip(*self.para_raios.values())
            fig.add_trace(go.Scatter(x=x_pr, y=y_pr, mode='markers+text', 
                                     name='Para-raios', text=list(self.para_raios.keys()),
                                     marker=dict(size=10, color='gray', symbol='x')))

        fig.update_layout(title="Geometria da Torre", xaxis_title="x (m)", yaxis_title="y (m)",
                          template="plotly_dark", showlegend=True, height=600)
        return fig
