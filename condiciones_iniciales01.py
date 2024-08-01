import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parametros para el espacio de partículas
particle_spacing = 0.012167  # Espaciado entre partículas del fluido
delta = particle_spacing  # Pequeño delta para no sobreponerse con las paredes

# Parámetros de la caja y del bloque de agua
box_length = 0.6
box_width0 = 0.2
box_width = round(box_width0 / particle_spacing) * particle_spacing
print('box_width modificate=',box_width)
box_height = 0.4

# Parámetros para el bloque de agua
water_length = 0.3
water_width = box_width - delta/2.  # pequeño delta para no sobreponerse con las paredes
water_height = 0.3
wall_particle_spacing = particle_spacing / 2  # Espaciado entre partículas de las paredes

# Constantes físicas
g = 9.82
gamma = 7.0
beta0 = 1.0  # Ajuste según sea necesario
beta = beta0 * 1.0
ht = water_height
c = beta * np.sqrt(2. * g * ht)
rho0 = 1000.
b = rho0 * c * c / gamma

# Crear partículas para el bloque de agua
x_water = np.arange(delta, water_length, particle_spacing)
y_water = np.arange(delta, water_width, particle_spacing)
z_water = np.arange(delta, water_height, particle_spacing)
xx_water, yy_water, zz_water = np.meshgrid(x_water, y_water, z_water)
water_particles = np.vstack([xx_water.ravel(), yy_water.ravel(), zz_water.ravel()]).T

# Crear partículas para las paredes de la caja
x_box = np.arange(0, box_length+wall_particle_spacing, wall_particle_spacing)
y0_box = np.arange(0, box_width+wall_particle_spacing, wall_particle_spacing)
y_box = np.arange(wall_particle_spacing, box_width, wall_particle_spacing)
z_box = np.arange(wall_particle_spacing, box_height+wall_particle_spacing, wall_particle_spacing)

wall_particles = []
# Paredes laterales
for y in [0, box_width]:
    for x in x_box:
        for z in z_box:
            wall_particles.append([x, y, z])
# Paredes frontal y trasera
for x in [0, box_length]:
    for y in y_box:
        for z in z_box:
            wall_particles.append([x, y, z])
# Piso
for x in x_box:
    for y in y0_box:
        wall_particles.append([x, y, 0])

wall_particles = np.array(wall_particles)

# Combinar partículas de la caja y del agua
num_fluid_particles = water_particles.shape[0]
num_box_particles = wall_particles.shape[0]
num_particles = num_fluid_particles + num_box_particles

# Inicializar IDs y velocidades
ids_fluid = np.arange(1, num_fluid_particles + 1)
ids_wall = np.arange(num_fluid_particles + 1, num_particles + 1)
velocities = np.zeros((num_particles, 3))

# Densidad, presión y masa de las partículas
#rho = np.full(num_particles, rho0)
dx = particle_spacing
rho = np.zeros(num_particles)
pressure = np.zeros(num_particles)
mass = np.zeros(num_particles)
for i in range(num_fluid_particles):
    rho[i] = rho0 * (1 + (rho0 * g * (ht - abs(water_particles[i, 2])) / b ))**(1.0 / gamma)
    pressure[i] = b * ((rho[i] / rho0)**gamma - 1.0)
    mass[i] = dx *dx *dx * rho[i]
    
# Para las partículas de pared, mantener rho0
rho[num_fluid_particles:] = rho0
pressure[num_fluid_particles:] = 0.0
mass[num_fluid_particles:] = dx * dx * dx * rho0

#mass = np.full(num_particles, dx * dx * dx * rho0)
#mass = np.full(num_particles, (4/3) * np.pi * dx*dx*dx * rho0)

# Presión y energía interna
internal_energy = np.full(num_particles, 357.1)

# Tipo de partícula y hsml
itype = np.full(num_particles, -1)
itype[:num_fluid_particles] = 1  # Partículas del fluido
hsml = np.full(num_particles, dx)

# Guardar los datos en un archivo
with open('snapshot_000', 'w') as f:
    f.write(f"0   0.0   {num_particles}   {num_fluid_particles}   {num_box_particles}\n")
    
    # Escribir partículas del fluido
    for i in range(num_fluid_particles):
        f.write(f"{ids_fluid[i]}   {water_particles[i, 0]:.6f}   {water_particles[i, 1]:.6f}   {water_particles[i, 2]:.6f}   "
                f"{velocities[i, 0]:.6f}   {velocities[i, 1]:.6f}   {velocities[i, 2]:.6f}   {mass[i]:.6f}   {rho[i]:.6f}   "
                f"{pressure[i]:.6f}   {internal_energy[i]:.6f}   {itype[i]}   {hsml[i]:.6f}\n")
    
    # Escribir partículas de las paredes
    for i in range(num_box_particles):
        j = num_fluid_particles + i
        f.write(f"{ids_wall[i]}   {wall_particles[i, 0]:.6f}   {wall_particles[i, 1]:.6f}   {wall_particles[i, 2]:.6f}   "
                f"{velocities[j, 0]:.6f}   {velocities[j, 1]:.6f}   {velocities[j, 2]:.6f}   {mass[j]:.6f}   {rho[j]:.6f}   "
                f"{pressure[j]:.6f}   {internal_energy[j]:.6f}   {itype[j]}   {hsml[j]:.6f}\n")

# Graficar las partículas en 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Paredes de la caja (puntos más pequeños y transparentes)
ax.scatter(wall_particles[:, 0], wall_particles[:, 1], wall_particles[:, 2], c='r', s=10, alpha=0.01, label='Box')


# Bloque de agua (puntos más grandes)
ax.scatter(water_particles[:, 0], water_particles[:, 1], water_particles[:, 2], c='b', s=50, label='Water')



# Definir los límites de los ejes para que tengan la misma escala
max_range = np.array([box_length, box_width, box_height]).max() / 2.0
mid_x = (x_box.max() + x_box.min()) * 0.5
mid_y = (y_box.max() + y_box.min()) * 0.5
mid_z = (z_box.max() + z_box.min()) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)

ax.set_box_aspect([1,1,1])  # Para asegurarse de que los ejes tengan la misma escala



ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
plt.show()
