from Dfit import *
from decimal import Decimal
# optional
tmin = 10 # Minimum timestep used for MSD calculation
tmax = 100 # Maximum timestep used for MSD calculation
m = 20 # Number of MSD values in fitting
nseg = 1000 # Number of segments
dt = 5 # dt*tunit is the timestep of the input trajectory
#path="/media/eniac/mdd1/paper_membranas/analisis/coef_diff_agua/"
#path="/media/eniac/mdd1/paper_membranas/snx_water_coef_difu/dcd/"
filename='memtox2dprotein.dat'
fz=filename

fout='mem_tox_coef'
res_3D = Dfit.Dcov(fz=fz,tmin=tmin,tmax=tmax,m=m,nseg=nseg,dt=dt,imgfmt='png')
res_3D.run_Dfit()




# For analysis:
tc = 5 # tc is connected to the timestep [in ps] by tc*dt. If tmin=1 and dt=1, then tc=7 corresponds to a 7ps timestep.

# If desired, for finite-size correction (Yeh, Hummer, JPCB 2004):
eta = 9.0e-4 # Viscosity in Pa*s
L = 7.7 # Box edge length in nm (cubic box!)
T = 300 # Temperature in Kelvin

res_3D.analysis(tc=tc)
#res_3D.finite_size_correction(T=T,eta=eta,L=L,tc=tc)

#0.004592890975343372 nm^2/ps
def nm_ps_to_cm2_s(value= 1.9689570536523056e-05):
    value=float( '%.2E' % Decimal(value))
    print (value)
    ps_to_seg=1e-12
    nm2_cm2=1e-14
    cm2=value*nm2_cm2
    print (cm2)
    cm2_s=cm2/ps_to_seg
    cm2_s='%.2E' % Decimal(cm2_s)
    print (cm2_s)
nm_ps_to_cm2_s()    