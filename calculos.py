import math
from decimal import Decimal
# nanometros_a_cm=1e-14
# slope=0.0003805936391483073
# DcdFreq=5000
# femto_a_seg=1e-15 
# dt=2
# dim=2
# D= (slope * nanometros_a_cm)/(2*dim *(dt*DcdFreq)*femto_a_seg)
# print (D) 
# print (1.90E-07/2)
mol=6.02e23
L=1/1000 
anstrom_a_cm=1e-08
anstrom3_a_cm3=1e-24
R=15* anstrom_a_cm #a

#radius 50.0 snx -- mem
D= 9.51E-08 #cm2/s
vol=108 * 102 * 14 * anstrom3_a_cm3
print (vol)
P=1/vol

ka=(4*math.pi)*D*R*P
print('%E' % Decimal(str(ka)))

