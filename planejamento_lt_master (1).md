# ⚡ LT-Master Analysis - Planejamento Completo v2.0

> Software de análise de linhas de transmissão baseado em Fuchs, Stevenson e rigor do ATPDraw

---

## 🎯 Visão Geral do Sistema

**Objetivo:** Sistema completo para análise elétrica e mecânica de linhas de transmissão de alta tensão, incluindo:
- Cálculo de parâmetros elétricos com efeito solo
- Modelagem de propagação e quadripolos
- Análise de fluxo de potência e estabilidade
- Compensação de reativos
- Limites térmicos e de operação
- Acoplamento entre circuitos

**Tecnologias Base:**
- Python 3.9+ (Backend/Engine)
- Streamlit (Interface responsiva)
- NumPy/SciPy (Álgebra matricial e otimização)
- Plotly/Matplotlib (Visualização interativa)
- Pandas (Manipulação de dados)

---

## 📦 Arquitetura de Módulos

### **FASE 1: Fundação Física e Geometria**

#### 🔷 Módulo 1: Dados e Estrutura da Torre
**Arquivo:** `engine/geometria.py`

**Componentes:**

1. **Banco de Dados de Condutores** (`data/condutores_db.json`)
   ```json
   {
     "Grosbeak": {
       "tipo": "CAA 636 MCM",
       "diametro_mm": 31.77,
       "rmg_mm": 12.65,
       "rdc_25C": 0.05328,
       "rac_25C": 0.05856,
       "rac_50C": 0.06197,
       "secao_mm2": 322.3,
       "corrente_max_A": 730,
       "peso_kg_km": 1087
     }
   }
   ```
   
   **Condutores inclusos:** Grosbeak, Cardinal, Rail, Pheasant, Dove, Bluejay, Drake, Tern, EHS 3/8", EHS 1/2"

2. **Classe GeometriaTorre**
   ```python
   class GeometriaTorre:
       def __init__(self, tipo_circuito='simples'):
           # tipo: 'simples' ou 'duplo'
           self.fases = {}      # {A, B, C: (x, y)}
           self.para_raios = {} # {PR1, PR2: (x, y)}
           self.circuito2 = {}  # Para circuito duplo
           
       def calcular_rmg_equivalente(self, num_condutores, espacamento):
           # Bundle: RMG_eq = ⁿ√(RMG × r^(n-1))
           # n = número de condutores
           # r = raio do círculo do bundle
           
       def calcular_dmd(self):
           # Matriz de distâncias mútuas D_ij
           # Para circuito duplo: incluir acoplamento
           
       def validar_geometria(self):
           # Verificar distâncias mínimas
           # Checagem de isolamento fase-terra
           
       def plotar_torre_2d(self):
           # Visualização esquemática
   ```

3. **Configuração de Bundle**
   - Suporte para 1 a 6 condutores por fase
   - Arranjos: horizontal, vertical, triangular, quadrado
   - Cálculo automático do raio equivalente do feixe

---

#### 🔷 Módulo 2: Parâmetros Elétricos da Linha
**Arquivo:** `engine/parametros.py`

##### **2.1 - Impedância de Série (Z)**

**Classe Principal:**
```python
class ImpedanciaLinha:
    def __init__(self, geometria, condutor, freq=60, rho_solo=100):
        self.geometria = geometria
        self.condutor = condutor
        self.freq = freq
        self.rho_solo = rho_solo  # Ω⋅m (NOVO - CRÍTICO)
```

**Funções:**

1. **Resistência com Correção Térmica e Efeito Pelicular**
   ```python
   def calcular_resistencia(self, temp_C):
       # Correção de temperatura
       alpha = 0.00403  # 1/°C para alumínio
       R_temp = self.R_dc * (1 + alpha * (temp_C - 25))
       
       # Efeito pelicular (simplificado)
       k_skin = self.R_ac / self.R_dc
       R_final = R_temp * k_skin
       
       return R_final  # Ω/km
   ```

2. **Matriz de Indutância Primitiva com Efeito Solo (Carson/Dubanton)**
   ```python
   def matriz_indutancia_primitiva(self):
       """
       Método de Dubanton para retorno pela terra
       """
       omega = 2 * pi * self.freq
       mu_0 = 4e-7 * pi  # H/m
       
       # Profundidade equivalente de retorno (Dubanton)
       D_e = 658.5 * sqrt(self.rho_solo / self.freq)  # metros
       
       n = len(condutores_total)  # fases + para-raios
       L_prim = np.zeros((n, n), dtype=complex)
       
       for i in range(n):
           for j in range(n):
               if i == j:
                   # Autoindutância
                   L_prim[i,i] = 2e-7 * ln(D_e / RMG_i)
               else:
                   # Indutância mútua
                   D_ij = distancia(i, j)
                   D_ij_img = distancia_imagem(i, j, D_e)
                   L_prim[i,j] = 2e-7 * ln(D_ij_img / D_ij)
       
       # Conversão para impedância
       Z_prim = R_prim + 1j * omega * L_prim
       return Z_prim  # Ω/km
   ```

3. **Redução de Kron (Eliminar Para-raios)**
   ```python
   def reducao_kron(self, Z_primitiva):
       """
       Reduzir matriz (n×n) para (3×3) eliminando para-raios
       Z_abc = Z_ff - Z_fg @ inv(Z_gg) @ Z_gf
       """
       n_fases = 3  # ou 6 para circuito duplo
       
       Z_ff = Z_primitiva[:n_fases, :n_fases]
       Z_fg = Z_primitiva[:n_fases, n_fases:]
       Z_gf = Z_primitiva[n_fases:, :n_fases]
       Z_gg = Z_primitiva[n_fases:, n_fases:]
       
       Z_abc = Z_ff - Z_fg @ np.linalg.inv(Z_gg) @ Z_gf
       return Z_abc
   ```

4. **Transposição (Opcional)**
   ```python
   def aplicar_transposicao(self, Z_abc):
       """
       Linha transposta: Z equilibrada
       Z_s = (Z_aa + Z_bb + Z_cc) / 3
       Z_m = (Z_ab + Z_bc + Z_ca) / 3
       """
       Z_proprio = (Z_abc[0,0] + Z_abc[1,1] + Z_abc[2,2]) / 3
       Z_mutuo = (Z_abc[0,1] + Z_abc[1,2] + Z_abc[2,0]) / 3
       
       Z_eq = np.array([
           [Z_proprio, Z_mutuo, Z_mutuo],
           [Z_mutuo, Z_proprio, Z_mutuo],
           [Z_mutuo, Z_mutuo, Z_proprio]
       ])
       return Z_eq
   ```

##### **2.2 - Admitância em Derivação (Y)**

**Funções:**

1. **Matriz de Potencial de Maxwell (com Método das Imagens)**
   ```python
   def matriz_potencial(self):
       """
       Capacitância usando método das imagens de Kelvin
       """
       epsilon_0 = 8.854e-12  # F/m
       
       n = len(condutores)
       P = np.zeros((n, n))
       
       for i in range(n):
           h_i = altura_condutor[i]
           r_i = raio_condutor[i]
           
           for j in range(n):
               if i == j:
                   # Potencial próprio
                   P[i,i] = (1/(2*pi*epsilon_0)) * ln(2*h_i / r_i)
               else:
                   # Potencial mútuo
                   D_ij = distancia(i, j)
                   D_ij_img = distancia_com_imagem(i, j)
                   P[i,j] = (1/(2*pi*epsilon_0)) * ln(D_ij_img / D_ij)
       
       return P  # m/F
   ```

2. **Matriz de Capacitância**
   ```python
   def matriz_capacitancia(self, P):
       """
       C = inv(P)
       Aplicar Kron para eliminar para-raios
       """
       C_primitiva = np.linalg.inv(P)
       C_abc = self.reducao_kron(C_primitiva)
       
       # Conversão para admitância (S/km)
       omega = 2 * pi * self.freq
       B_abc = omega * C_abc * 1e-3  # S/km
       
       return B_abc
   ```

3. **Condutância (Perdas)**
   ```python
   def calcular_condutancia(self):
       """
       G ≈ ω × C × tan(δ)
       tan(δ) ≈ 1e-6 (perdas dielétricas + corona)
       """
       tan_delta = 1e-6
       G = self.freq * 2*pi * self.C_abc * tan_delta
       return G  # S/km
   ```

##### **2.3 - Componentes Simétricas**

**Arquivo:** `engine/sequencias.py`

```python
def transformacao_fortescue(Z_abc, Y_abc):
    """
    Conversão ABC → 012
    """
    a = np.exp(1j * 2*pi/3)  # Operador de Fortescue
    
    A = (1/3) * np.array([
        [1,  1,   1  ],
        [1,  a**2, a ],
        [1,  a,   a**2]
    ])
    
    A_inv = np.array([
        [1,  1,   1  ],
        [1,  a,   a**2],
        [1,  a**2, a ]
    ])
    
    Z_012 = A @ Z_abc @ A_inv
    Y_012 = A @ Y_abc @ A_inv
    
    return {
        'Z0': Z_012[0,0],
        'Z1': Z_012[1,1],
        'Z2': Z_012[2,2],
        'Y0': Y_012[0,0],
        'Y1': Y_012[1,1],
        'Y2': Y_012[2,2]
    }
```

**⚠️ PONTO CRÍTICO:** A impedância de sequência zero (Z0) depende **fortemente** da resistividade do solo. Sem este parâmetro, o cálculo de faltas fase-terra será impreciso.

---

### **FASE 2: Propagação e Análise Operacional**

#### 🔷 Módulo 3: Parâmetros de Propagação
**Arquivo:** `engine/propagacao.py`

**Cálculos:**

1. **Impedância Característica**
   ```python
   def impedancia_caracteristica(Z, Y):
       """
       Zc = sqrt(Z / Y)
       """
       Zc_seq = {}
       for seq in ['0', '1', '2']:
           Zc_seq[seq] = np.sqrt(Z[seq] / Y[seq])
       return Zc_seq  # Ω
   ```

2. **Constante de Propagação**
   ```python
   def constante_propagacao(Z, Y):
       """
       γ = sqrt(Z × Y) = α + jβ
       α = atenuação (Np/km)
       β = constante de fase (rad/km)
       """
       gamma = np.sqrt(Z * Y)
       alpha = gamma.real
       beta = gamma.imag
       
       return {
           'gamma': gamma,
           'alpha': alpha,  # Np/km
           'beta': beta     # rad/km
       }
   ```

3. **Velocidade e Comprimento de Onda**
   ```python
   def parametros_onda(beta, freq):
       omega = 2 * pi * freq
       v = omega / beta  # km/s
       lambda_onda = 2*pi / beta  # km
       
       # Converter para unidades práticas
       v_pct_luz = (v / 300000) * 100  # % da velocidade da luz
       
       return {
           'velocidade_km_s': v,
           'velocidade_pct_luz': v_pct_luz,
           'comprimento_onda_km': lambda_onda
       }
   ```

4. **Modelo de Quadripolo (ABCD)**
   ```python
   def matriz_abcd(Z, Y, comprimento_km, tipo='longa'):
       """
       Calcula ABCD conforme o modelo da linha
       """
       l = comprimento_km
       
       if tipo == 'curta' or l < 80:
           # Linha Curta
           A = D = 1
           B = Z * l
           C = 0
           
       elif tipo == 'media' or l < 250:
           # Linha Média (modelo π)
           A = D = 1 + (Z*l * Y*l) / 2
           B = Z * l
           C = Y*l * (1 + (Z*l * Y*l) / 4)
           
       else:
           # Linha Longa (modelo exato com funções hiperbólicas)
           gamma = np.sqrt(Z * Y)
           Zc = np.sqrt(Z / Y)
           
           A = D = np.cosh(gamma * l)
           B = Zc * np.sinh(gamma * l)
           C = np.sinh(gamma * l) / Zc
       
       return np.array([[A, B], [C, D]])
   ```

---

#### 🔷 Módulo 4: Fluxo de Potência e Estabilidade
**Arquivo:** `engine/fluxo_potencia.py`

**Análises:**

1. **Cálculo de Terminal Receptor**
   ```python
   def calcular_receptor(Vs, delta_s, P_recebida, FP, ABCD):
       """
       Dado Vs, Ps, calcula Vr e Ir
       Método iterativo (Newton-Raphson simplificado)
       """
       A, B, C, D = ABCD[0,0], ABCD[0,1], ABCD[1,0], ABCD[1,1]
       
       # Estimativa inicial
       Vr = Vs / abs(A)
       
       for iter in range(10):  # Máx 10 iterações
           Sr = P_recebida / FP + 1j*P_recebida*tan(acos(FP))
           Ir = conj(Sr / Vr)
           
           # Atualizar Vr pela equação do quadripolo
           Vr_novo = (Vs - B*Ir) / A
           
           if abs(Vr_novo - Vr) < 0.01:
               break
           Vr = Vr_novo
       
       Is = C*Vs + D*Ir
       
       return {
           'Vr': Vr,
           'Ir': Ir,
           'Is': Is,
           'delta_r': angle(Vr)
       }
   ```

2. **Regulação de Tensão**
   ```python
   def regulacao_tensao(Vr_vazio, Vr_carga):
       """
       Reg% = ((|Vr_vazio| - |Vr_carga|) / |Vr_carga|) × 100
       """
       reg = ((abs(Vr_vazio) - abs(Vr_carga)) / abs(Vr_carga)) * 100
       return reg
   ```

3. **Eficiência da Transmissão**
   ```python
   def eficiencia(Ps, Pr):
       """
       η = (Pr / Ps) × 100
       """
       return (Pr.real / Ps.real) * 100
   ```

4. **Curva P×δ (Estabilidade Estática)**
   ```python
   def curva_potencia_angulo(Vs, Vr, Z_eq):
       """
       P = (|Vs|⋅|Vr| / |Z|) × sin(δ)
       """
       delta_range = np.linspace(-90, 90, 181)  # graus
       P = []
       
       for delta in delta_range:
           delta_rad = np.radians(delta)
           P_delta = (abs(Vs) * abs(Vr) / abs(Z_eq)) * np.sin(delta_rad)
           P.append(P_delta)
       
       P_max = abs(Vs) * abs(Vr) / abs(Z_eq)  # Limite estático (δ = 90°)
       
       return {
           'delta': delta_range,
           'potencia': P,
           'P_max': P_max
       }
   ```

5. **Potência Natural (SIL)**
   ```python
   def calcular_sil(V_nominal_kV, Zc):
       """
       SIL = |V|² / |Zc|  (MW)
       """
       V_ln = V_nominal_kV / sqrt(3)  # kV linha-neutro
       SIL_MW = (V_ln**2 / abs(Zc)) * 1000  # MW por fase
       SIL_total = 3 * SIL_MW
       
       return SIL_total
   ```

---

### **FASE 3: Compensação e Limites Operacionais**

#### 🔷 Módulo 5: Compensação de Reativos
**Arquivo:** `engine/compensacao.py`

**Tipos de Compensação:**

1. **Compensação Shunt (Derivação)**
   ```python
   class CompensacaoShunt:
       def dimensionar_reator(self, V_nominal, efeito_ferranti_pc):
           """
           Reator para limitar sobretensão em vazio
           Q_reator = V² / X_reator
           """
           V = V_nominal
           Q_necessaria = (V**2 / abs(self.Zc)) * (efeito_ferranti_pc / 100)
           X_reator = V**2 / Q_necessaria
           
           return {
               'Q_Mvar': Q_necessaria,
               'X_ohm': X_reator,
               'localizacao': 'terminal_receptor'
           }
       
       def dimensionar_capacitor(self, V_nominal, Q_compensar):
           """
           Capacitor para suporte de tensão em carga
           """
           X_capacitor = V_nominal**2 / Q_compensar
           
           return {
               'Q_Mvar': Q_compensar,
               'X_ohm': X_capacitor
           }
   ```

2. **Compensação Série**
   ```python
   class CompensacaoSerie:
       def calcular_grau_compensacao(self, X_linha, k_comp):
           """
           X_compensado = X_linha × (1 - k_comp)
           k_comp típico: 0.3 a 0.7
           """
           X_capacitor_serie = -k_comp * X_linha
           X_resultante = X_linha + X_capacitor_serie
           
           return {
               'k_compensacao': k_comp,
               'X_capacitor': X_capacitor_serie,
               'X_linha_compensada': X_resultante,
               'aumento_capacidade_pc': (1/(1-k_comp) - 1) * 100
           }
   ```

3. **Análise Comparativa**
   ```python
   def comparar_com_sem_compensacao(self, caso_base, caso_compensado):
       """
       Plota perfis de tensão, perdas e capacidade
       """
       melhorias = {
           'reducao_perdas_pc': ...,
           'melhoria_regulacao_pc': ...,
           'aumento_capacidade_MW': ...
       }
       return melhorias
   ```

---

#### 🔷 Módulo 6: Ampacidade e Limites Térmicos (NOVO - CRÍTICO)
**Arquivo:** `premium/ampacidade.py`

**Baseado na IEEE 738:**

```python
class AmpacidadeLinha:
    def __init__(self, condutor, temp_max_C=75):
        self.condutor = condutor
        self.temp_max = temp_max_C
    
    def calcular_limite_termico(self, temp_amb, vento_m_s, radiacao_W_m2):
        """
        Balanço térmico:
        q_joule = q_conveccao + q_radiacao - q_solar
        
        I² × R = h_conv × A × ΔT + σ × ε × A × (T⁴ - T_amb⁴) - α_solar × A × S
        """
        # Resistência na temperatura máxima
        R_max = self.condutor.calcular_R(self.temp_max)
        
        # Área superficial do condutor
        D = self.condutor.diametro_mm / 1000  # metros
        A = pi * D  # m²/m
        
        # Coeficiente de convecção (função do vento)
        h_conv = self.calcular_conveccao(vento_m_s, D)
        
        # Perdas por radiação
        sigma = 5.67e-8  # W/(m²⋅K⁴)
        emissividade = 0.5  # Alumínio oxidado
        q_rad = sigma * emissividade * A * (
            (self.temp_max + 273.15)**4 - (temp_amb + 273.15)**4
        )
        
        # Ganho solar
        absortividade = 0.5
        q_solar = absortividade * A * radiacao_W_m2
        
        # Perdas por convecção
        q_conv = h_conv * A * (self.temp_max - temp_amb)
        
        # Corrente máxima
        q_total_dissipavel = q_conv + q_rad - q_solar
        I_max = sqrt(q_total_dissipavel / R_max)
        
        return {
            'I_max_A': I_max,
            'P_max_MW': sqrt(3) * self.V_nominal * I_max / 1e6,
            'condicoes': {
                'temp_amb': temp_amb,
                'vento': vento_m_s,
                'radiacao': radiacao_W_m2
            }
        }
    
    def calcular_conveccao(self, vento, diametro):
        """
        Correlação de Nusselt para cilindro em fluxo cruzado
        """
        # Simplificado: h ≈ 10 + 6×sqrt(vento)
        h = 10 + 6 * np.sqrt(vento)
        return h  # W/(m²⋅K)
```

**Gráfico de Capabilidade:**
```python
def plotar_curva_capabilidade(self, comprimentos_km):
    """
    Mostra 3 limites:
    1. Térmico (horizontal)
    2. Queda de tensão (Vr > 0.95 pu)
    3. Estabilidade (δ < 30°)
    """
    P_termico = self.I_max * sqrt(3) * V_nominal / 1e6
    
    P_queda_tensao = []
    P_estabilidade = []
    
    for L in comprimentos_km:
        # Calcular limite por queda de tensão
        # Calcular limite por estabilidade
        ...
    
    plt.plot(comprimentos_km, [P_termico]*len(comprimentos_km), label='Térmico')
    plt.plot(comprimentos_km, P_queda_tensao, label='Regulação 5%')
    plt.plot(comprimentos_km, P_estabilidade, label='Estabilidade')
    plt.xlabel('Comprimento (km)')
    plt.ylabel('Potência Máxima (MW)')
    plt.legend()
```

---

### **FASE 4: Módulos Avançados (Premium)**

#### 🔷 Módulo 7: Acoplamento entre Linhas (NOVO - IMPORTANTE)
**Arquivo:** `premium/acoplamento.py`

**Para Circuitos Duplos ou Linhas Paralelas:**

```python
class AcoplamentoMutuo:
    def __init__(self, linha1, linha2, distancia_entre_eixos):
        self.linha1 = linha1
        self.linha2 = linha2
        self.dist = distancia_entre_eixos
    
    def calcular_impedancia_mutua(self):
        """
        Z_mutua entre circuitos paralelos
        Importante para sequência zero
        """
        # Matriz ampliada (6×6 para 2 circuitos trifásicos)
        Z_primitiva_6x6 = self.montar_matriz_completa()
        
        # Redução de Kron mantendo as 6 fases
        Z_6x6 = self.reducao_kron(Z_primitiva_6x6)
        
        # Extrair blocos
        Z11 = Z_6x6[0:3, 0:3]  # Circuito 1
        Z22 = Z_6x6[3:6, 3:6]  # Circuito 2
        Z12 = Z_6x6[0:3, 3:6]  # Acoplamento
        
        return {
            'Z_circuito1': Z11,
            'Z_circuito2': Z22,
            'Z_mutua': Z12
        }
    
    def fator_desequilibrio_sequencia_zero(self):
        """
        Quanto a corrente de sequência zero de um circuito
        induz tensão no outro
        """
        Z_mutua_0 = self.transformar_para_sequencia_zero(self.Z12)
        return abs(Z_mutua_0)
```

---

#### 🔷 Módulo 8: Efeito Corona
**Arquivo:** `premium/corona.py`

```python
def tensao_critica_disruptiva(condutor, altitude_m, condicao_superficie):
    """
    Fórmula de Peek
    """
    r = condutor.raio_mm / 1000  # metros
    D = condutor.DMD  # metros (distância entre fases)
    
    # Fator de rugosidade
    m_0 = {
        'liso': 1.0,
        'oxidado': 0.98,
        'sujo': 0.87
    }[condicao_superficie]
    
    # Densidade relativa do ar
    delta = exp(-altitude_m / 8150)
    
    # Gradiente crítico (kV/cm)
    E_0 = 21.1 * m_0 * delta * (1 + 0.3 / sqrt(delta * r * 100))
    
    # Tensão crítica fase-neutro (kV)
    V_critica = E_0 * r * 100 * ln(D / r)
    
    # Tensão de linha
    V_critica_linha = V_critica * sqrt(3)
    
    return V_critica_linha
```

---

#### 🔷 Módulo 9: Campos Eletromagnéticos
**Arquivo:** `premium/campos_em.py`

```python
def campo_eletrico_solo(x, y, geometria, tensoes_fases):
    """
    Método das imagens para calcular E(x,y) no solo
    """
    E_total = 0 + 0j
    
    for fase, (x_fase, y_fase) in geometria.fases.items():
        V_fase = tensoes_fases[fase]
        
        # Contribuição direta
        r = sqrt((x - x_fase)**2 + (y - y_fase)**2)
        E_direto = V_fase / (2*pi*epsilon_0 * r)
        
        # Contribuição da imagem
        r_img = sqrt((x - x_fase)**2 + (y + y_fase)**2)
        E_imagem = -V_fase / (2*pi*epsilon_0 * r_img)
        
        E_total += E_direto + E_imagem
    
    return abs(E_total)  # kV/m

def campo_magnetico_solo(x, y, geometria, correntes_fases):
    """
    Lei de Biot-Savart
    """
    mu_0 = 4e-7 * pi
    B_total = 0 + 0j
    
    for fase, (x_fase, y_fase) in geometria.fases.items():
        I_fase = correntes_fases[fase]
        r = sqrt((x - x_fase)**2 + (y - y_fase)**2)
        
        B = (mu