import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parámetros para el cañón aleatorio
ncols, nrows = 100, 100  # Número de columnas y filas
cellsize = 1.  # Tamaño de cada celda en las coordenadas
xllcorner, yllcorner = 0.0, 0.0  # Coordenadas de la esquina inferior izquierda
angle_degrees_y = 10  # Ángulo de inclinación en grados en dirección Y
angle_degrees_x = 5   # Ángulo de inclinación en grados en dirección X

# Conversión de ángulo a radianes
angle_radians_y = np.radians(angle_degrees_y)
angle_radians_x = np.radians(angle_degrees_x)

# Generación del cañón aleatorio usando una combinación de funciones gaussianas
x = np.linspace(xllcorner, xllcorner + (ncols - 1) * cellsize, ncols)
y = np.linspace(yllcorner, yllcorner + (nrows - 1) * cellsize, nrows)
X, Y = np.meshgrid(x, y)

# Modificación para hacer los picos más anchos
ancho = 0.005  # Reducido de 0.01 a 0.005 para hacer los picos más anchos
xpico1 = 25    # Posición del pico 1
xpico2 = 75    # Posición del pico 2
alturaPicos = 50.  # Altura de los picos
profundidadValle = 20   # Profundidad del valle
xvalle = 50. # Posición del valle
Z = alturaPicos * (np.exp(-ancho * ((X - xllcorner - xpico1)**2)) + np.exp(-ancho * ((X - xllcorner - xpico2)**2)))  # Picos laterales más anchos
Z -= profundidadValle * np.exp(-ancho * ((X - xllcorner - xvalle)**2))  # Valle central menos pronunciado
Z += np.random.normal(0, 0.2, (nrows, ncols))  # Añadiendo ruido suave para realismo

# Añadir inclinaciones al terreno
# Inclinación en dirección Y
Z += Y * np.tan(angle_radians_y)
# Inclinación en dirección X (segunda inclinación)
Z += X * np.tan(angle_radians_x)

# Crear archivo de datos similar al archivo leído por el código proporcionado
filename = 'ascii_cañon_inclinado_diagonal.txt'
with open(filename, 'w') as file:
    file.write(f"ncols         {ncols}\n")
    file.write(f"nrows         {nrows}\n")
    file.write(f"xllcorner     {xllcorner}\n")
    file.write(f"yllcorner     {yllcorner}\n")
    file.write(f"cellsize      {cellsize}\n")
    file.write("NODATA_value  -9999\n")
    np.savetxt(file, Z, fmt='%.2f')

# Visualización en 3D con malla
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='k')

# Configuración de la gráfica
ax.set_title(f'Cañón Aleatorio en 3D con Inclinación Diagonal ({angle_degrees_y}° en Y, {angle_degrees_x}° en X)')
ax.set_xlabel('Coordenada X')
ax.set_ylabel('Coordenada Y')
ax.set_zlabel('Altitud')

# Configurar la relación de aspecto de los ejes para que sea la misma en X, Y y Z
ax.set_box_aspect([np.ptp(x), np.ptp(y), np.ptp(Z)])  # Relación de aspecto

plt.show()

