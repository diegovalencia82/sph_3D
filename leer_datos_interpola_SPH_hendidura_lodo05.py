import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

from collections import defaultdict

##### =================== FUNCIONES=======================

def interpolar_puntos1D(X,Y,Z, num_puntos=100, metodo='cubic'):

    # Interpolación para obtener una cuadrícula uniforme si es necesario
    #n_interp = 200  # Número de puntos de interpolación en cada eje

    x_interp = np.linspace(np.min(X), np.max(X), num_puntos)
    y_interp = np.linspace(np.min(Y), np.max(Y), num_puntos)
    X_interp, Y_interp = np.meshgrid(x_interp, y_interp)
    ddx = x_interp[3]-x_interp[2]
    ddy = y_interp[3]-y_interp[2]
    print('cellsize boundary  dx =',ddx,'cellsize boundary dy =',ddy)
    # Interpolar los datos Z en la cuadrícula uniforme
    Z_interp = griddata((X.flatten(), Y.flatten()), Z.flatten(), (X_interp, Y_interp), method='cubic')

    return X_interp, Y_interp, Z_interp, ddx


def interpolar_puntos(puntos, num_puntos=100, metodo='cubic'):
    """
    Interpola los puntos de la matriz en una cuadrícula regular.

    Parámetros:
    - hendidura: numpy array con las coordenadas (X, Y, Z) de los puntos de la hendidura.
    - num_puntos: número de puntos en cada eje de la cuadrícula (opcional, por defecto 100).
    - metodo: método de interpolación (opcional, por defecto 'cubic').

    Retorna:
    - Xh_grid: cuadrícula de coordenadas X.
    - Yh_grid: cuadrícula de coordenadas Y.
    - Zh_grid: valores interpolados de Z en la cuadrícula.
    """
    # Crear una cuadrícula para interpolar los puntos
    xh_grid = np.linspace(np.min(puntos[:, 0]), np.max(puntos[:, 0]), num_puntos)
    yh_grid = np.linspace(np.min(puntos[:, 1]), np.max(puntos[:, 1]), num_puntos)
    Xh_grid, Yh_grid = np.meshgrid(xh_grid, yh_grid)

    # Interpolar los valores Z de la hendidura sobre la cuadrícula
    Zh_grid = griddata((puntos[:, 0], puntos[:, 1]), puntos[:, 2], (Xh_grid, Yh_grid), method=metodo)
    
    return Xh_grid, Yh_grid, Zh_grid





def generar_cilindro_con_corte(x0, y0, z0, r, altura, delta, x1, y1, z1, x2, y2, z2, theta):
    # Calcular el vector de dirección del plano de corte a partir de los dos puntos
    dir_vector = np.array([x2 - x1, y2 - y1, z2 - z1])
    
    # Normalizar el vector de dirección
    norm_dir_vector = dir_vector / np.linalg.norm(dir_vector)

    # Crear un vector que representa la dirección en el plano xy usando el ángulo theta
    theta_rad = np.radians(theta)
    rotation_vector = np.array([np.cos(theta_rad), np.sin(theta_rad), 0])

    # Generar la cuadrícula de puntos del cilindro
    Xcilindro = []
    Ycilindro = []
    Zcilindro = []

    # Generar puntos en el cilindro desde z0 hasta z0 + altura
    for z in np.arange(z0, z0 + altura, delta):
        for t in np.arange(0, 2 * np.pi, delta / r):
            x = x0 + r * np.cos(t)
            y = y0 + r * np.sin(t)
            Xcilindro.append(x)
            Ycilindro.append(y)
            Zcilindro.append(z)

    # Cerrar la parte inferior del cilindro (base llena de puntos)
    for x in np.arange(x0 - r, x0 + r, delta):
        for y in np.arange(y0 - r, y0 + r, delta):
            if np.sqrt((x - x0)**2 + (y - y0)**2) <= r:
                Xcilindro.append(x)
                Ycilindro.append(y)
                Zcilindro.append(z0)

    # Aplicar el corte inclinado
    Xcortado = []
    Ycortado = []
    Zcortado = []
    X_eliminados = []  # Lista para los puntos eliminados
    Y_eliminados = []
    Z_eliminados = []

    # Recorrer los puntos del cilindro y comprobar si están por debajo del plano de corte
    for x, y, z in zip(Xcilindro, Ycilindro, Zcilindro):
        # Ecuación del plano: n • (P - P0) = 0
        # Donde n es el vector normal al plano
        # P es el punto (x, y, z)
        # P0 es uno de los puntos del plano, por ejemplo, (x1, y1, z1)

        # Calcular el vector normal al plano a partir de los puntos
        normal_vector = np.cross(norm_dir_vector, rotation_vector)
        normal_vector = normal_vector / np.linalg.norm(normal_vector)  # Normalizar el vector normal

        # Ecuación del plano: n • (P - P0) = 0
        d = np.dot(normal_vector, np.array([x - x1, y - y1, z - z1]))

        if d <= 0:  # Mantener los puntos por debajo del plano de corte
            Xcortado.append(x)
            Ycortado.append(y)
            Zcortado.append(z)
        else:  # Guardar los puntos eliminados
            X_eliminados.append(x)
            Y_eliminados.append(y)
            Z_eliminados.append(z)

    return (np.array(Xcortado), np.array(Ycortado), np.array(Zcortado),
            np.array(X_eliminados), np.array(Y_eliminados), np.array(Z_eliminados))


def obtener_paredes_cilindro(Xcortado, Ycortado, Zcortado, x0, y0, r, tolerancia=1e-3):
    """
    Filtra los puntos de las paredes del cilindro (sin incluir la base).
    
    Parámetros:
    - Xcortado, Ycortado, Zcortado: vectores con las coordenadas de los puntos del cilindro cortado.
    - x0, y0: centro del cilindro.
    - r: radio del cilindro.
    - tolerancia: margen de error para filtrar los puntos de las paredes.
    
    Retorna:
    - xc, yc, zc: vectores que representan los puntos en las paredes del cilindro.
    """
    
    # Inicializar listas para las coordenadas de las paredes del cilindro
    xc, yc, zc = [], [], []
    
    # Encontrar la altura mínima (base) del cilindro
    z_min = np.min(Zcortado)
    
    # Recorrer todos los puntos del cilindro cortado
    for x, y, z in zip(Xcortado, Ycortado, Zcortado):
        # Calcular la distancia radial al centro del cilindro
        distancia_al_centro = np.sqrt((x - x0)**2 + (y - y0)**2)
        
        # Filtrar los puntos que están en las paredes del cilindro (cerca del radio) y no en la base
        if np.abs(distancia_al_centro - r) <= tolerancia and z > z_min:
            xc.append(x)
            yc.append(y)
            zc.append(z)
    
    return np.array(xc), np.array(yc), np.array(zc)

def obtener_puntos_mas_altos_cilindro(Xcortado, Ycortado, Zcortado, x0, y0, r, tolerancia=1e-3):
    """
    Filtra los puntos más altos de cada columna del cilindro, sin incluir la base.
    
    Parámetros:
    - Xcortado, Ycortado, Zcortado: vectores con las coordenadas de los puntos del cilindro cortado.
    - x0, y0: centro del cilindro.
    - r: radio del cilindro.
    - tolerancia: margen de error para filtrar los puntos de las paredes.
    
    Retorna:
    - xc, yc, zc: vectores que representan los puntos más altos en las paredes del cilindro.
    """
    
    # Inicializar un diccionario para almacenar los puntos más altos en cada columna (x, y)
    puntos_mas_altos = {}
    
    # Encontrar la altura mínima (base) del cilindro
    z_min = np.min(Zcortado)
    
    # Recorrer todos los puntos del cilindro cortado
    for x, y, z in zip(Xcortado, Ycortado, Zcortado):
        # Calcular la distancia radial al centro del cilindro
        distancia_al_centro = np.sqrt((x - x0)**2 + (y - y0)**2)
        
        # Filtrar los puntos que están en las paredes del cilindro (cerca del radio) y no en la base
        if np.abs(distancia_al_centro - r) <= tolerancia and z > z_min:
            # Crear una clave única para cada columna (x, y) usando una tupla redondeada a la tolerancia
            clave_columna = (round(x, 5), round(y, 5))
            
            # Si la columna aún no está en el diccionario o el nuevo punto es más alto, actualizar el valor
            if clave_columna not in puntos_mas_altos or puntos_mas_altos[clave_columna][2] < z:
                puntos_mas_altos[clave_columna] = (x, y, z)
    
    # Extraer las coordenadas más altas en cada columna
    xc, yc, zc = [], [], []
    for x, y, z in puntos_mas_altos.values():
        xc.append(x)
        yc.append(y)
        zc.append(z)
    
    return np.array(xc), np.array(yc), np.array(zc)




def llenar_volumen_cilindro_completo(Xcortado, Ycortado, Zcortado, x0, y0, r, delta, tolerancia=1e-3):
    """
    Llena todo el volumen del cilindro cortado con partículas de fluido hasta la misma altura
    de los puntos más altos del corte.

    Parámetros:
    - Xcortado, Ycortado, Zcortado: vectores con las coordenadas de los puntos del cilindro cortado.
    - x0, y0: centro del cilindro.
    - r: radio del cilindro.
    - delta: distancia entre las partículas de fluido.
    - tolerancia: margen de error para filtrar los puntos de las paredes.

    Retorna:
    - xfluido, yfluido, zfluido: vectores con las coordenadas de las partículas de fluido.
    """

    # Obtener los puntos más altos de cada columna en el cilindro cortado
    xc_pared, yc_pared, zc_pared = obtener_puntos_mas_altos_cilindro(Xcortado, Ycortado, Zcortado, x0, y0, r, tolerancia)

    # Crear una malla uniforme para recorrer todas las posiciones dentro del radio del cilindro
    x_values = np.arange(x0 - r, x0 + r, delta)
    y_values = np.arange(y0 - r, y0 + r, delta)

    # Inicializar listas para las coordenadas del fluido
    xfluido, yfluido, zfluido = [], [], []

    # Recorrer todos los puntos dentro del radio del cilindro
    for x in x_values:
        for y in y_values:
            # Verificar que el punto esté dentro del radio del cilindro
            if np.sqrt((x - x0)**2 + (y - y0)**2) <= r-delta/2.:
                # Encontrar la altura máxima permitida en esta columna (x, y)
                # Calcular las distancias entre (x, y) y los puntos de las paredes
                diferencias_x = np.abs(xc_pared - x)  # Diferencias en x
                indices_cercanos = np.where(diferencias_x == np.min(diferencias_x))  # Índices de los puntos más cercanos en x
            
                # Obtener las alturas correspondientes
                z_max = zc_pared[indices_cercanos]  # Alturas de los puntos más cercanos en x
                
                #Tomar la máxima altura entre los más cercanos en x
                z_max = np.max(z_max) if z_max.size > 0 else z_min  # Asegúrate de manejar el caso donde no hay puntos
                
                # Generar partículas de fluido desde la base hasta la altura z_max
                for z in np.arange(np.min(Zcortado)+delta*1.5, z_max, delta):
                    xfluido.append(x)
                    yfluido.append(y)
                    zfluido.append(z)

    return np.array(xfluido), np.array(yfluido), np.array(zfluido)








def rotar_punto(x0, y0, xp, yp, beta):
    # Convertir el ángulo beta a radianes
    #beta_rad = np.radians(beta)
    
    # Trasladar el punto para que (x0, y0) sea el origen
    x_trans = xp - x0
    y_trans = yp - y0
    
    # Aplicar la rotación
    x_rotado = x_trans * np.cos(beta) - y_trans * np.sin(beta)
    y_rotado = x_trans * np.sin(beta) + y_trans * np.cos(beta)
    
    # Trasladar de vuelta el punto a la posición original
    x_final = x_rotado + x0
    y_final = y_rotado + y0
    
    # Retornar las coordenadas rotadas
    return x_final, y_final
    


############ ================ FIN DE FUNCIONES ================================


# Leer los datos del archivo ASCII
filename = 'ascii_cañon_inclinado_diagonal.txt'
with open(filename, 'r') as file:
    header = {}
    for i in range(6):
        line = file.readline().strip().split()
        header[line[0]] = float(line[1])
    data = np.loadtxt(file)

ncols = int(header['ncols'])
nrows = int(header['nrows'])
xllcorner = header['xllcorner']
yllcorner = header['yllcorner']
cellsize = header['cellsize']

print('cellsize=',cellsize)

# Crear malla de coordenadas X, Y
x = np.linspace(xllcorner, xllcorner + (ncols - 1) * cellsize, ncols)
y = np.linspace(yllcorner, yllcorner + (nrows - 1) * cellsize, nrows)
XA, YA = np.meshgrid(x, y)

# Crear una copia de los datos originales para modificarlos
ZA = np.copy(data)

X, Y, Z, ddx = interpolar_puntos1D(XA, YA, ZA, num_puntos=200, metodo='cubic')


n_filas, n_columnas = X.shape

print("Forma de X_interp:", X.shape)
print("Forma de Y_interp:", Y.shape)
print("Forma de Z_interp:", Z.shape)


# =================== Posición de la hendidura ===============================

x0, y0 = 60, 75  # Coordenadas del centro de la hendidura
r =5            # Radio del trozo de círculo

# Variables para almacenar el punto más bajo y alto antes de la hendidura
min_altitud_original = np.inf
max_altitud_original = -np.inf
min_pos_original = (None, None)
max_pos_original = (None, None)

# Identificar el punto más bajo dentro del radio antes de aplicar la hendidura
for i in range(n_filas):
    for j in range(n_columnas):
        # Calcular la distancia desde el punto actual hasta el centro de la hendidura
        distance = np.sqrt((X[i, j] - x0)**2 + (Y[i, j] - y0)**2)
        if distance < r:
            # Verificar si este punto es el más bajo en la superficie original
            if Z[i, j] < min_altitud_original:
                min_altitud_original = Z[i, j]
                min_pos_original = (X[i, j], Y[i, j])
            # Verificar si este punto es el más bajo en la superficie original
            if Z[i, j] > max_altitud_original:
                max_altitud_original = Z[i, j]
                max_pos_original = (X[i, j], Y[i, j])

# Imprimir el punto más bajo en la superficie original
print(f"El punto más bajo en la superficie original está en las coordenadas X={min_pos_original[0]}, Y={min_pos_original[1]} con una altitud de Z={min_altitud_original}")
X0min = min_pos_original[0]
Y0min = min_pos_original[1]
Z0min = min_altitud_original
# Imprimir el punto más alto en la superficie original
print(f"El punto más alto en la superficie original está en las coordenadas X={max_pos_original[0]}, Y={max_pos_original[1]} con una altitud de Z={max_altitud_original}")
X0max = max_pos_original[0]
Y0max = max_pos_original[1]
Z0max = max_altitud_original

ht = max_altitud_original-min_altitud_original
print('Diferencia ht =', ht)

# ================================================================================================

# Inicializar arrays vacíos para la hendidura con la misma forma que X, Y, Z
X_hendidura = np.zeros((n_filas, n_columnas))
Y_hendidura = np.zeros((n_filas, n_columnas))
Z_hendidura = np.zeros((n_filas, n_columnas))

p_corte = []

X_lodo = []
Y_lodo = []
Z_lodo = []

delta_z = cellsize * np.sqrt(2)  # Espaciado entre los puntos de lodo

for i in range(n_filas):
    for j in range(n_columnas):
        X_hendidura[i, j] = X[i, j]
        Y_hendidura[i, j] = Y[i, j]
        Z_hendidura[i, j] = Z[i, j]
        
        d1 = np.sqrt((X[i, j] - x0)**2 + (Y[i, j] - y0)**2)
        if d1 < r:
            if r-ddx < d1 < r+ddx:
                p_corte.append([ X[i, j],Y[i, j], Z[i, j] ])
            if d1 == 0.0:
                X_hendidura[i, j] = 0.0#X0max
                Y_hendidura[i, j] = 0.0#Y0max
            else:
                d2 = r - d1
                k = d2 / d1
                X_hendidura[i, j] = -99999.#X[i, j] - k * (x0 - X[i, j])
                Y_hendidura[i, j] = -99999.#Y[i, j] - k * (y0 - Y[i, j])           
                v1x = X0max - x0
                v1y = Y0max - y0
                v2x = X_hendidura[i, j] - x0
                v2y = Y_hendidura[i, j] - y0
                cross_product = v1x * v2y - v1y * v2x
                dot_product = v1x * v2x + v1y * v2y
                theta = np.arctan2(cross_product, dot_product)

                if np.abs(theta) > np.pi/2 :
                    if theta > 0:
                        beta = -( 2*np.abs(theta) - np.pi )
                    else:
                        beta = ( 2*np.abs(theta) - np.pi )
                        
                    #X_hendidura[i, j], Y_hendidura[i, j] = rotar_punto(x0, y0, X_hendidura[i, j], Y_hendidura[i, j], beta)
                

topox = []
topoy = []
topoz = []
k = -1
for i in range(n_filas):
    for j in range(n_columnas):
        if X_hendidura[i,j] != -99999. and Y_hendidura[i,j] != -99999.:
            x = X_hendidura[i,j]
            y = Y_hendidura[i,j]
            z = Z_hendidura[i,j]
            topox.append(x)
            topoy.append(y)
            topoz.append(z)
            
topox = np.array(topox)
topoy = np.array(topoy)
topoz = np.array(topoz)


# Parámetros del cilindro
#x0, y0, z0 = 0, 0, 0  # Centro y altura inicial del cilindro
r = 5  # Radio del cilindro
altura = ht  # Altura del cilindro
delta = 0.5  # Espaciado entre puntos

# Puntos del plano de corte

x1, y1, z1 = X0min, Y0min, Z0min  # Primer punto del plano
x2, y2, z2 = X0max, Y0max, Z0max  # Segundo punto del plano

# Ángulo de dirección en grados
theta = 90  # Dirección del plano de corte en el plano xy

# Generar el cilindro y aplicar el corte
Xcortado, Ycortado, Zcortado, X_eliminados, Y_eliminados, Z_eliminados = generar_cilindro_con_corte(x0, y0, Z0min, r, altura, ddx, x1, y1, z1, x2, y2, z2, theta)

# Unir la falla Xcortado  con la topografia topox
X_topo = np.concatenate((Xcortado, topox))
Y_topo = np.concatenate((Ycortado, topoy))
Z_topo = np.concatenate((Zcortado, topoz))


# Llenar el cilindro cortado con fluido
deltaf = 2*ddx
# Generar las partículas de fluido en el cilindro cortado
X_fluido, Y_fluido, Z_fluido = llenar_volumen_cilindro_completo(Xcortado, Ycortado, Zcortado, x0, y0, r, deltaf, deltaf/5.)

# ======================== # Asignar propiedades a las partículas del fluido========


# Constantes físicas
g = 9.82
gamma = 7.0
beta0 = 1.0  # Ajuste según sea necesario
beta = beta0 * 1.0
ht = max_altitud_original #-min_altitud_original
c = beta * np.sqrt(2. * g * ht)
rho0 = 1000.0  # Densidad base del fluido
b = rho0 * c * c / gamma
print('Constantes Físicas')
print('g = ',g)
print('gamma =', 7.0)
print('beta = beta0 * 1.0 = ',beta)
print('ht = Y0max - Y0c = ',ht)
print('c = beta * np.sqrt(2. * g * ht) = ',c)
print('rho0 =',rho0)
print('b = rho0 * c * c / gamma = ',b)

# IDs y tipos de partículas
num_wall_particles = len(X_topo)
num_fluid_particles = len(X_fluido)
num_particles = num_wall_particles + num_fluid_particles
ids_fluid = np.arange(1, num_fluid_particles + 1)  # IDs para partículas del fluido
ids_wall = np.arange(num_fluid_particles + 1, num_particles + 1)  # IDs para partículas de la topografía
print('==========================Variables para el SPH ')
print('num_wall_particles=', num_wall_particles, 'num_fluid_particles=', num_fluid_particles, 'num_particles=', num_particles,len(X_topo))
print('delta wall = ',ddx, 'delta fluid = ',deltaf)
print('==========================')

mass_wall = rho0 * deltaf**3        # Masa para las partículas de la topografía
mass_fluid =rho0 * ddx**3        # Masa para las partículas del fluido

# Inicializar arrays para las propiedades de las partículas
velocities = np.zeros((num_particles, 3))  # Inicializa la velocidad para todas las partículas
rho = np.zeros(num_particles)
pressure = np.zeros(num_particles)
mass = np.zeros(num_particles)
internal_energy = np.full(num_particles, 357.1)
itype = np.full(num_particles, -1)
hsml = np.full(num_particles, deltaf)

itype[:num_fluid_particles] = 1  # Tipo 1 para las partículas del fluido
for i in range(num_fluid_particles):
    rho[i] = rho0 * (1 + (rho0 * g * (ht - abs(Z_fluido[i])) / b ))**(1.0 / gamma)
    pressure[i] = b * ((rho[i] / rho0)**gamma - 1.0)
    mass[i] = rho[i] * deltaf**3


# Asignar propiedades a las partículas de la topografía
rho[num_fluid_particles:] = rho0
pressure[num_fluid_particles:] = 0.0
mass[num_fluid_particles:] = mass_wall
itype[num_fluid_particles:] = -1  # Tipo 0 para las partículas de la topografía

# Guardar los datos en un archivo
with open('snapshot_000', 'w') as f:
    f.write(f"0   0.0   {num_particles}   {num_fluid_particles}   {num_wall_particles}\n")

    # Escribir partículas del fluido
    for i in range(num_fluid_particles):
        f.write(f"{ids_fluid[i]}   {X_fluido[i]:.6f}   {Y_fluido[i]:.6f}   {Z_fluido[i]:.6f}   "
                f"{velocities[i, 0]:.6f}   {velocities[i, 1]:.6f}   {velocities[i, 2]:.6f}   {mass[i]:.6f}   {rho[i]:.6f}   "
                f"{pressure[i]:.6f}   {internal_energy[i]:.6f}   {itype[i]}   {hsml[i]:.6f}\n")
    
    # Escribir partículas de la topografía   
    for i in range(num_wall_particles):
        j = num_fluid_particles + i
        f.write(f"{ids_wall[i]}   {X_topo[i]:.6f}   {Y_topo[i]:.6f}   {Z_topo[i]:.6f}   "
                f"{velocities[j, 0]:.6f}   {velocities[j, 1]:.6f}   {velocities[j, 2]:.6f}   {mass[j]:.6f}   {rho[j]:.6f}   "
                f"{pressure[j]:.6f}   {internal_energy[j]:.6f}   {itype[j]}   {hsml[j]:.6f}\n")


# =======================================   GRAFICAS   =============================


# OPCIÓN 1: Graficar superficie de la montaña con puntos verdes y lodo con puntos cafés
fig1 = plt.figure(figsize=(10, 7))
ax1 = fig1.add_subplot(111, projection='3d')

# Graficar los puntos
#ax1.scatter(X_hendidura, Y_hendidura, Z_hendidura, c='g', marker='o', s=0.1, label='Puntos Superficie')
ax1.scatter(topox, topoy, topoz, c='g', marker='o', s=0.1, label='Puntos Superficie')
ax1.scatter(Xcortado, Ycortado, Zcortado, c='r', marker='o', s=0.5, label='Puntos Superficie')
#ax1.scatter(xc, yc, zc, c='b', marker='o', s=5, label='Puntos Fluido')
ax1.scatter(X_fluido, Y_fluido, Z_fluido, c='b', marker='o', s=5, label='Puntos Fluido')
#ax1.scatter(X_eliminados, Y_eliminados, Z_eliminados, c='r', marker='o', s=0.1, label='Puntos Superficie')
ax1.scatter(X0max, Y0max, Z0max, c='r', marker='o', s=5, label='Punto máximo')
ax1.scatter(X0min, Y0min, Z0min, c='r', marker='o', s=5, label='Punto mínimo')
ax1.set_title('Superficie de la Montaña y Puntos de Lodo')
ax1.set_xlabel('Coordenada X')
ax1.set_ylabel('Coordenada Y')
ax1.set_zlabel('Altitud (Z)')
ax1.legend()
plt.show()
