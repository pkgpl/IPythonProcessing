# CDP-sorted data
input_filename="marm_gain_decon_agc_sort.sd"

# NMO-corrected data
output_filename="marm_gain_decon_agc_sort_nmo.sd"

vmin=1500.	# minimum velocity
dv=100		# velocity increment
nv=41		# number of velocity samples (max. vel = vmin + (nv-1)*dv )

cmp_start=25	# first cdp number
cmp_step=50	# cdp step for velocity analysis
ncmp=11		# number of cdps to pick velocity

max_stretch=10.0	# max. allowable nmo stretch in %

R = 4 # t0: time step (samples)
L = 8 # semblance window half length (samples)

