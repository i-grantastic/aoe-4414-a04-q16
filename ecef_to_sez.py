# ecef_to_sez.py
#
# Usage: python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km
#   converts ECEF coordinates to SEZ coordinates

# Parameters:
#   o_x_km: x-component of the origin in km
#   o_y_km: y-component of the origin in km
#   o_z_km: z-component of the origin in km
#   x_km: x-component of the position in km
#   y_km: y-component of the position in km
#   z_km: z-component of the position in km

# Output:
#   SEZ coordinates
#
# Written by Grant Chapman
# Other contributors: None

# import Python modules
import math # math module
import sys # argv
import numpy as np

# constants
R_E_KM = 6378.137
E_E    = 0.081819221456

## calculate demoninator
# (eccentricity, latitude in radians)
def calc_denom(ecc, lat_rad):
  return math.sqrt(1.0-ecc**2.0 * math.sin(lat_rad)**2.0)

# initialize script arguments
o_x_km = float('nan')
o_y_km = float('nan')
o_z_km = float('nan')
x_km = float('nan')
y_km = float('nan')
z_km = float('nan')

# parse script arguments
if len(sys.argv) == 7:
  o_x_km = float(sys.argv[1])
  o_y_km = float(sys.argv[2])
  o_z_km = float(sys.argv[3])
  x_km = float(sys.argv[4])
  y_km = float(sys.argv[5])
  z_km = float(sys.argv[6])
else:
  print(\
    'Usage: '\
    'python3 ecef_to_sez.py o_x_km o_y_km o_z_km x_km y_km z_km'\
  )
  exit()

### script below this line ###

# calculate longitude
lon_rad = math.atan2(o_y_km, o_x_km)

# initialize lat_rad, r_lon_km, r_z_km
lat_rad = math.asin(o_z_km/math.sqrt(o_x_km**2+o_y_km**2+o_z_km**2))
r_lon_km = math.sqrt(o_x_km**2+o_y_km**2)
prev_lat_rad = float('nan')

# iteratively find latitude
c_E = float('nan')
count = 0
while (math.isnan(prev_lat_rad) or abs(lat_rad-prev_lat_rad)>10e-7) and count<5:
  denom = calc_denom(E_E,lat_rad)
  c_E = R_E_KM/denom
  prev_lat_rad = lat_rad
  lat_rad = math.atan((o_z_km+c_E*(E_E**2)*math.sin(lat_rad))/r_lon_km)
  count = count+1

# calculate hae
hae_km = r_lon_km/math.cos(lat_rad)-c_E

# rotation matrices
r_i_z = [[math.sin(lat_rad), 0, -math.cos(lat_rad)],
         [0, 1, 0],
         [math.cos(lat_rad), 0, math.sin(lat_rad)]]
r_i_y = [[math.cos(lon_rad), math.sin(lon_rad), 0],
         [-math.sin(lon_rad), math.cos(lon_rad), 0],
         [0, 0, 1]]

# ECEF vector
ecef_vector = [x_km-o_x_km, y_km-o_y_km, z_km-o_z_km]

# first rotation
first_rotation = np.matmul(r_i_y, ecef_vector)

# second rotation
second_rotation = np.matmul(r_i_z, first_rotation)

# SEZ vector
s_km = second_rotation[0]
e_km = second_rotation[1]
z_km = second_rotation[2]

# print
print(s_km)
print(e_km)
print(z_km)