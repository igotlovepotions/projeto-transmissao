import streamlit as st
import numpy as np
import pandas as pd
import json
import plotly.graph_objects as go
from engine.geometria import GeometriaTorre
from engine.parametros import ParametrosLinha
from engine.sequencias import transformacao_fortescue
from engine.propagacao import PropagacaoLinha
from engine.fluxo_potencia import FluxoPotencia
from premium.ampacidade import AmpacidadeLinha

# Configuração da Página
st.set_page_config(page_title="LT-Master Analysis v2.0", layout="wide", initial_sidebar_state="expanded")

# CSS customizado para visual premium
st.markdown("""
    <style>
    .main {
        background-color: #0e1117;
    }
    .stMetric {
        background-color: #1e2130;
        padding: 15px;
        border-radius: 10px;
        border: 1px solid #4a4a4a;
    }
    h1, h2, h3 {
        color: #00d4ff;
    }
    </style>
    """, unsafe_allow_html=True)

st.title("⚡ LT-Master Analysis - Dashboard")
st.markdown("---")

# Carregar banco de dados de condutores
with open('data/condutores_db.json', 'r') as f:
    condutores_db = json.load(f)

# Sidebar - Configurações Gerais
st.sidebar.header("🛠️ Configurações do Projeto")
v_nom_ll = st.sidebar.number_input("Tensão Nominal (kV LL)", value=230.0, step=10.0)
freq = st.sidebar.selectbox("Frequência (Hz)", [60, 50], index=0)
rho_solo = st.sidebar.number_input("Resistividade do Solo (Ω.m)", value=100.0)
comprimento = st.sidebar.number_input("Comprimento da Linha (km)", value=100.0)

# Layout da Tower Geometry
st.sidebar.header("📐 Geometria da Torre")
tipo_torre = st.sidebar.selectbox("Tipo de Torre", ["Simples", "Duplo"])
cond_fase_nome = st.sidebar.selectbox("Condutor de Fase", list(condutores_db.keys()))
num_cond = st.sidebar.slider("Nº condutores por feixe", 1, 6, 2)
espac_feixe = st.sidebar.number_input("Espaçamento do feixe (m)", value=0.45)

cond_pr_nome = st.sidebar.selectbox("Cabo Para-raio", ["EHS 3/8\"", "EHS 1/2\""])

# Criar Geometria
geo = GeometriaTorre(tipo_circuito=tipo_torre.lower())
if tipo_torre == "Simples":
    geo.adicionar_fase('A', -7.0, 20.0)
    geo.adicionar_fase('B', 0.0, 20.0)
    geo.adicionar_fase('C', 7.0, 20.0)
    geo.adicionar_para_raio('PR1', -4.0, 28.0)
    geo.adicionar_para_raio('PR2', 4.0, 28.0)
else:
    # Exemplo circuito duplo
    geo.adicionar_fase('A1', -6.0, 18.0)
    geo.adicionar_fase('B1', -7.5, 23.0)
    geo.adicionar_fase('C1', -6.0, 28.0)
    geo.adicionar_fase('A2', 6.0, 18.0)
    geo.adicionar_fase('B2', 7.5, 23.0)
    geo.adicionar_fase('C2', 6.0, 28.0)
    geo.adicionar_para_raio('PR1', 0.0, 35.0)

# Cálculos
cond_fase = condutores_db[cond_fase_nome]
cond_pr = condutores_db[cond_pr_nome]

params = ParametrosLinha(geo, cond_fase, cond_pr, num_cond, espac_feixe, freq, rho_solo)
z_prim = params.matriz_impedancia_primitiva()
z_abc = params.reducao_kron(z_prim)
y_abc = params.calcular_admitancia()

seq = transformacao_fortescue(z_abc, y_abc)

prop = PropagacaoLinha(
    {'0': seq['Z0'], '1': seq['Z1'], '2': seq['Z2']},
    {'0': seq['Y0'], '1': seq['Y1'], '2': seq['Y2']},
    freq
)
zc = prop.calcular_impedancia_caracteristica()
gamma = prop.calcular_constante_propagacao()
abcd = prop.matriz_abcd('1', comprimento)

# Dashboard Layout
col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("📍 Geometria da Torre")
    st.plotly_chart(geo.plotar_torre_2d(), use_container_width=True)

with col2:
    st.subheader("🔢 Parâmetros de Sequência (Ω/km & S/km)")
    data_seq = {
        "Parâmetro": ["Resistência (R)", "Reatância (X)", "Susceptância (B)"],
        "Sequência 1": [seq['Z1'].real, seq['Z1'].imag, seq['Y1'].imag],
        "Sequência 0": [seq['Z0'].real, seq['Z0'].imag, seq['Y0'].imag]
    }
    st.table(pd.DataFrame(data_seq))
    
    st.subheader("🌊 Propagação (Seq 1)")
    onda = prop.calcular_parametros_onda(gamma['1']['beta'])
    st.metric("Impedância Característica (Zc)", f"{abs(zc['1']):.2f} Ω")
    st.metric("Comprimento de Onda", f"{onda['comprimento_onda_km']:.1f} km")
    st.metric("Velocidade de Propagação", f"{onda['velocidade_pct_luz']:.1f} % da luz")

st.markdown("---")
st.subheader("🔌 Análise de Operação e Carga")

col_op1, col_op2 = st.columns(2)

with col_op1:
    p_carga = st.slider("Potência de Carga (MW)", 0.0, 1000.0, 200.0)
    fp = st.slider("Fator de Potência", 0.7, 1.0, 0.95)
    
    fluxo = FluxoPotencia(abcd, v_nom_ll)
    res_fluxo = fluxo.calcular_terminal_receptor(1.0, p_carga, fp)
    sil = fluxo.calcular_sil(zc['1'])

    st.metric("Tensão no Receptor (Vr)", f"{res_fluxo['Vr_pu']:.3f} pu")
    st.metric("Regulação de Tensão", f"{res_fluxo['Regulacao']:.2f} %")
    st.metric("Eficiência da Transmissão", f"{res_fluxo['Eficiencia']:.2f} %")
    st.metric("Potência Natural (SIL)", f"{sil:.2f} MW")

with col_op2:
    st.subheader("⚖️ Ampacidade e Limites Térmicos")
    t_amb = st.slider("Temperatura Ambiente (°C)", 0, 45, 25)
    vento = st.slider("Velocidade do Vento (m/s)", 0.1, 10.0, 0.6)
    
    # Estimativa de R na temp max
    r_max = cond_fase['rac_50C'] * 1.15 # Approx p/ 75C
    termico = AmpacidadeLinha(cond_fase['diametro_mm'], r_max)
    res_termico = termico.calcular_limite_termico(t_amb, vento)
    
    st.metric("Corrente Máxima (I_max)", f"{res_termico['I_max_A']:.1f} A")
    i_op = abs(res_fluxo['Is_a'])
    pct_termico = (i_op / res_termico['I_max_A']) * 100
    st.progress(min(pct_termico/100, 1.0))
    st.write(f"Utilização Térmica: {pct_termico:.1f}%")

st.markdown("---")
st.info("💡 Projeto desenvolvido com base no planejamento LT-Master v2.0")
