# CDP-sorted data
input_filename="trc_gain_bpf_sd_agc_sort.sd"

# NMO-corrected data
output_filename="trc_gain_bpf_sd_agc_sort_nmo.sd"

vmin=5000.	# minimum velocity
dv=200		# velocity increment
nv=51		# number of velocity samples (max. vel = vmin + (nv-1)*dv )

cmp_start=205	# first cdp number
cmp_step=5	# cdp step for velocity analysis
ncmp=13		# number of cdps to pick velocity

max_stretch=10.0	# max. allowable nmo stretch in %

R = 4 # t0: time step (samples)
L = 8 # semblance window half length (samples)

