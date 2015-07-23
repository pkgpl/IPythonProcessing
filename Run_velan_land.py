## CDP-sorted data
#input_filename="marm_gain_decon_agc_sort.sd"
#
## NMO-corrected data
#output_filename="marm_gain_decon_agc_sort_nmo_test.sd"
#
#vmin=1500.	# minimum velocity
#dv=100		# velocity increment
#nv=41		# number of velocity samples (max. vel = vmin + (nv-1)*dv )
#
#cmp_start=25	# first cdp number
#cmp_step=100	# cdp step for velocity analysis
#ncmp=6		# number of cdps to pick velocity
#
#max_stretch=10.0	# max. allowable nmo stretch in %
#
#R = 4 # t0: time step (samples)
#L = 8 # semblance window half length (samples)

from Par_velan_land import *

# do not edit codes below
cmp_end=cmp_start+(ncmp-1)*cmp_step	# last cdp number

print("cmp_start: %s"%cmp_start)
print("cmp_step: %s"%cmp_step)
print("cmp_end: %s"%cmp_end)
print("ncmp: %s"%ncmp)
print("CAUTION: At least two picks per panel are required")

# read CMP-sorted data
import pkprocess.pkbase as pk
su = pk.read(input_filename)

# check cmp range
cmpu=pk.get_key_unique(su,'cdp')
cmpmin=cmpu.min()
cmpmax=cmpu.max()
if cmp_start < cmpmin:
	print("cmp_start too small! (min.cmp=%s)"%cmpmin)
	import sys
	sys.exit(1)
if cmp_end > cmpmax:
	print("cmp_end too large! (max.cmp=%s)"%cmpmax)
	import sys
	sys.exit(1)

# Velocity analysis and NMO
from pkprocess.pkvelan import vel_picking_nmo
so = vel_picking_nmo(su,vmin,dv,nv,cmp_start,cmp_end,cmp_step,max_stretch,R,L)

# write NMO-corrected data
print("saving file: %s"%output_filename)
so.write(output_filename)
