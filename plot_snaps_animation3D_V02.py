import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.gridspec import GridSpec

# Variable de selección para mostrar o no las partículas con itype = -1
sel = 0  # 1 para mostrar las partículas con itype = -1, 0 para no mostrarlas

def read_snapshot(filepath):
    # Leer datos del archivo snapshot
    with open(filepath, 'r') as file:
        lines = file.readlines()
        # Extraer el tiempo del segundo dato de la primera línea que contiene metadatos
        time = float(lines[0].split()[1])
        # Saltar la primera línea que contiene metadatos
        lines = lines[1:]
        # Extraer las posiciones x, y, z, las velocidades vx, vy, vz y otros datos de cada partícula de cada línea
        data = np.array([[float(val) for val in line.split()] for line in lines])
    return time, data

def calculate_velocity(data):
    # Calcular la magnitud de la velocidad para cada partícula en 3D
    velocity = np.sqrt(data[:, 4]**2 + data[:, 5]**2 + data[:, 6]**2)  # Considerando las columnas 4, 5 y 6 como las velocidades vx, vy, vz
    return velocity

def update(frame):
    ax_3d.clear()
    ax_xz.clear()
    
    # Limites de la gráfica 3D
    ax_3d.set_xlim(-0.1, 0.8)  # Ajustar límites x según tus necesidades
    ax_3d.set_ylim(-0.1, 1.1)  # Ajustar límites y según tus necesidades
    ax_3d.set_zlim(-0.1, 0.8)  # Ajustar límites z según tus necesidades

    # Establecer los nombres de los ejes para la gráfica 3D
    ax_3d.set_xlabel('X')
    ax_3d.set_ylabel('Y')
    ax_3d.set_zlabel('Z')
    
    # Limites de la gráfica 2D (x-z)
    ax_xz.set_xlim(-0.1, 0.8)
    ax_xz.set_ylim(-0.1, 0.8)

    # Establecer los nombres de los ejes para la gráfica 2D (x-z)
    ax_xz.set_xlabel('X')
    ax_xz.set_ylabel('Z')
    
    # Filtrar partículas si sel es 0
    if sel == 0:
        data = frame[1][frame[1][:, 11] != -1]
    else:
        data = frame[1]

    # Calcular la velocidad de las partículas en el frame actual
    velocity = calculate_velocity(data)
    
    # Crear una escala de colores que vaya de azul a rojo
    norm = Normalize(vmin=0, vmax=np.sqrt(2*9.8*0.4))
    cmap = plt.get_cmap('coolwarm')
    sm = ScalarMappable(norm=norm, cmap=cmap)
    
    # Obtener el color de cada partícula según su velocidad
    colors = sm.to_rgba(velocity)
    
    # Dibujar las partículas en el frame actual con el color correspondiente en la gráfica 3D
    ax_3d.scatter(data[:, 1], data[:, 2], data[:, 3], color=colors,s=25,alpha=0.5)

    # Dibujar las partículas en el frame actual con el color correspondiente en la gráfica 2D (x-z)
    ax_xz.scatter(data[:, 1], data[:, 3], color=colors,s=25,alpha=0.5)

    if sel == 1:
        # Seleccionar las partículas con itype = -1 y pintarlas de color gris en ambas gráficas
        gray_particles = frame[1][frame[1][:, 11] == -1]
        ax_3d.scatter(gray_particles[:, 1], gray_particles[:, 2], gray_particles[:, 3], color='gray', s=10,alpha=0.01)
        ax_xz.scatter(gray_particles[:, 1], gray_particles[:, 3], color='gray',s=10,alpha=0.01)
    
    # Mostrar el tiempo en la gráfica 3D
    ax_3d.text2D(0.05, 0.95, f't = {frame[0]:.2f}', transform=ax_3d.transAxes, ha='left', va='top', fontsize=12)

# Definir la carpeta donde se encuentran los archivos
folder = './'

# Crear una figura con GridSpec para organizar las subtramas
fig = plt.figure(figsize=(12, 6))
gs = GridSpec(1, 2, width_ratios=[2, 1])

# Crear una subtrama para la animación en 3D
ax_3d = fig.add_subplot(gs[0], projection='3d')

# Crear una subtrama para la vista 2D (x-z)
ax_xz = fig.add_subplot(gs[1])

# Lista para almacenar los datos de todas las instantáneas
data_frames = []

# Iterar sobre los nombres de archivo del tipo 'snapshot_0001', 'snapshot_0002', ..., 'snapshot_2000
for i in range(1,500):
    # Formatear el nombre de archivo con el número de snapshot
    filename = f'snapshot_{i:04d}'
    # Crear la ruta completa al archivo
    filepath = os.path.join(folder, filename)
    # Verificar si el archivo existe
    if os.path.exists(filepath):
        # Leer datos del archivo snapshot y almacenarlos en la lista
        time, data = read_snapshot(filepath)
        data_frames.append((time, data))
    else:
        print(f'El archivo {filename} no existe en la carpeta especificada.')

# Crear la animación
anim = FuncAnimation(fig, update, frames=data_frames, interval=75)  # Intervalo de 5 milisegundos entre cada frame

# Guardar la animación como un archivo .mp4
anim.save('animacion_3D.mp4', writer='ffmpeg')

plt.show()
