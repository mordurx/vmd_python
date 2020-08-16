from trajectory import Trajectory

sim1_dcd= "/mnt/e/rep0_tox_mem500.dcd"
sim1_psf=" rep0.mem.tox.psf"
traj=Trajectory(sim1_dcd,sim1_psf)
print (traj)