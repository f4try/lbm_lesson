WIDTH = 256
L_ch = 0.5E-2
DX = L_ch / (WIDTH - 1)
VMU_H2 = 2.0984E-5
VMU_AIR = 4.2179E-5
FTAO = 1.
VMU_LAT = (FTAO - 0.5) / 3.
SCALE = VMU_AIR / VMU_LAT
DT = pow(DX, 2.) / SCALE
print(DT)