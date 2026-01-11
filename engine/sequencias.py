import numpy as np
from typing import Dict

def transformacao_fortescue(z_abc: np.ndarray, y_abc: np.ndarray) -> Dict[str, complex]:
    """
    Realiza a transformação de Fortescue para obter as componentes de sequência 012.
    Para sistemas equilibrados (transpostos), as componentes de sequência são os
    elementos da diagonal da matriz transformada.
    """
    a = np.exp(1j * 2 * np.pi / 3)  # Operador de Fortescue
    
    # Matriz de transformação T (A)
    # Vabc = T * V012
    # T = [[1, 1, 1], [1, a^2, a], [1, a, a^2]]
    t = np.array([
        [1, 1, 1],
        [1, a**2, a],
        [1, a, a**2]
    ])
    
    t_inv = np.linalg.inv(t)
    
    # Z012 = T_inv * Zabc * T
    z_012 = t_inv @ z_abc @ t
    y_012 = t_inv @ y_abc @ t
    
    # Extraindo as sequências (diagonal)
    # 0 = Zero, 1 = Positiva, 2 = Negativa
    return {
        'Z0': z_012[0, 0],
        'Z1': z_012[1, 1],
        'Z2': z_012[2, 2],
        'Y0': y_012[0, 0],
        'Y1': y_012[1, 1],
        'Y2': y_012[2, 2]
    }

def aplicar_transposicao(matriz: np.ndarray) -> np.ndarray:
    """
    Aplica a simplificação de linha transposta:
    Z_s = média da diagonal
    Z_m = média dos elementos fora da diagonal
    """
    n = len(matriz)
    diag = np.mean(np.diag(matriz))
    
    # Soma de todos os elementos menos a diagonal
    off_diag_sum = np.sum(matriz) - np.sum(np.diag(matriz))
    off_diag = off_diag_sum / (n**2 - n)
    
    mat_eq = np.full((n, n), off_diag, dtype=complex)
    np.fill_diagonal(mat_eq, diag)
    
    return mat_eq
