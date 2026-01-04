import matplotlib.pyplot as plt
import math

# Definindo as distâncias dos cabos
cables = {
    'A': {'altura': 28, 'distancia': -5},
    'B': {'altura': 34, 'distancia': 0},
    'C': {'altura': 28, 'distancia': 5},
    'PR1': {'altura': 39, 'distancia': 0},
   # 'PR2': {'altura': 39, 'distancia': 7},
}

# Criar o gráfico
fig, ax = plt.subplots(figsize=(10, 6))

# Desenhar os cabos com uma única cor (preto) e manter as anotações
for cabo, dados in cables.items():
    x = dados['distancia']
    y = dados['altura']
    ax.plot(x, y, 'o',
            markersize=12, label=f"Cabo {cabo}", color='black') # Definido cor preta
    # Adicionar anotações com o nome do cabo
    ax.annotate(cabo,
                xy=(x, y),
                xytext=(5, 5), textcoords='offset points',
                fontsize=10, fontweight='bold')

    # Adicionar linhas de projeção para os eixos (esquerda e para baixo) - visibilidade melhorada
    ax.plot([x, x], [0, y], 'k--', linewidth=1.0, alpha=0.8) # Linha vertical para o eixo X
    ax.plot([0, x], [y, y], 'k--', linewidth=1.0, alpha=0.8) # Linha horizontal para o eixo Y


# Definir limites dos eixos baseados nos dados reais
# Adicionar margem de 20% para melhor visualização
alturas = [d['altura'] for d in cables.values()]

margem_y = (max(alturas) - min(alturas)) * 0.2

# Eixo Y começando em 0 e com margem superior
ax.set_ylim(0, max(alturas) + margem_y)

# Eixo X centralizado em 0, com limites de -20 a 20 e ticks de 2 em 2
ax.set_xlim(-20, 20)
ax.set_xticks(range(-20, 21, 2))

# Adicionar detalhes visuais
ax.set_xlabel('Posição Horizontal (m)', fontsize=12)
ax.set_ylabel('Distância Vertical (m)', fontsize=12) # Alterado o rótulo do eixo Y
ax.set_title('Visualização de Posição dos Cabos', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3, linestyle='--')
ax.legend(loc='best')

# Adicionar linha de referência no zero
ax.axvline(x=0, color='black', linestyle='-', linewidth=0.5, alpha=0.5)
ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5, alpha=0.5)

plt.tight_layout()
plt.show()

# Calcular e imprimir as distâncias entre os cabos (usando Teorema de Pitágoras)
print("\nDistâncias euclidianas entre os cabos:")
points = {name: (data['distancia'], data['altura']) for name, data in cables.items()}

# Create a list of cable names to iterate through for pairs
cable_names = list(points.keys())

for i in range(len(cable_names)):
    for j in range(i + 1, len(cable_names)):
        name1 = cable_names[i]
        name2 = cable_names[j]
        p1 = points[name1]
        p2 = points[name2]

        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        distance = math.hypot(dx, dy) # sqrt(dx**2 + dy**2)
        print(f"Distância entre {name1} e {name2}: {distance:.2f} m")