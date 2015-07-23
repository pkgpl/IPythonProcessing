import numpy as np
import matplotlib.pyplot as plt
import pkbase as pk

def perc_clip_val(data,perc=100):
	# return clipping value by percent clip, does not apply clipping
	# data: array to clip
	# perc: percent value
	# output: min/max clip values
	mperc=(100.-perc)*0.01*0.5
	tmp=np.sort(data.flatten())
	minloc=int(len(tmp)*mperc)
	maxloc=len(tmp)-minloc-1
	clipmin=tmp[minloc]
	clipmax=tmp[maxloc]
	return clipmin,clipmax

def perc_clip(data,perc=100):
	# clip data
	# data: array to clip
	# perc=100: percent value, clip min/max (100-perc)/2 percent
	# output: clipped array
	if perc == 100:
		return data
	clipmin,clipmax=perc_clip_val(data,perc)
	return data.clip(clipmin,clipmax)

def subplot_xylabel(i,xlabel,ylabel):
	fs="large"
	plt.xlabel(xlabel,fontsize=fs)
	if i==0:
		plt.ylabel(ylabel,fontsize=fs)
	else:
		plt.gca().set_yticklabels([])

def plot_comp(su_tuple,plot="image",figsize=[0,0],perc=100,cmap="gray_r",fill=True,scale=1,dx=0,key=False):
	# Compare SeismicTrace plots
	# su_tuple: tuple of traces to plot
	# plot="image": plot type, image ("image") or wiggle ("wiggle")
	# figsize: matplotlib figure size
	# perc=100: percent clip
	# cmap="gray_r": matplotlib colormap
	# fill=True: fill values greater than zero (wiggle)
	# scale=1: scale wiggle trace plots (wiggle)
	# dx=0: receiver interval (specfk)
	# key: x-axis label (image)
	if figsize==[0,0]:
		figsize=[10,10]
		if len(su_tuple) > 2:
			figsize=[15,10]
	plt.figure(figsize=figsize)
	ndata=len(su_tuple)
	for i,su in enumerate(su_tuple):
		plt.subplot(1,ndata,i+1)
		if plot=="image":
			plot_image(su,perc=perc,cmap=cmap,subplot=True,key=key)
			if key:
				subplot_xylabel(i,key.upper(),"Time (s)")
			else:
				subplot_xylabel(i,"Trace number","Time (s)")
		elif plot=="wiggle":
			plot_wiggle(su,perc=perc,fill=fill,scale=scale,subplot=True)
			subplot_xylabel(i,"Trace number","Time (s)")
		elif plot=="specfx":
			specfx(su,perc=perc,cmap=cmap,subplot=True)
			subplot_xylabel(i,"Trace number","Frequency (Hz)")
		elif plot=="specfk":
			specfk(su,perc=perc,cmap=cmap,dx=dx,subplot=True)
			subplot_xylabel(i,"Wavenumber","Frequency (Hz)")

		plt.gca().xaxis.tick_top()
		plt.gca().xaxis.set_label_position("top")
	plt.tight_layout()

def plot_image(self,figsize=[5,10],perc=100,cmap="gray_r",ratio=0,f2=0,d2=1,subplot=False,key=False):
	# plot image traces
	# figsize=[5,10]: matplotlib figure size [inch]
	# perc=100: percent clip
	# cmap="gray_r": matplotlib colormap
	# ratio=0: figure x,y aspect ratio
	plotdata=self.data.T
	plotdata=perc_clip(plotdata,perc)
	print("min=%s max=%s"%(plotdata.min(),plotdata.max()))

	ntr=pk.get_ntr(self)
	ns=pk.get_ns(self)
	dt=pk.get_dt(self)

	xmin=f2
	xmax=f2+(ntr-1)*d2

	if key in pk.SU_KEY_LIST:
		xmin=self.header[key][0]
		xmax=self.header[key][-1]

	tmax=(ns-1)*dt
	if ratio==0:
		ratio="auto"
	else:
		ratio=xmax/tmax*ratio

	if not subplot: plt.figure(figsize=figsize)
	plt.imshow(plotdata,aspect=ratio,extent=[xmin,xmax,tmax,0],cmap=cmap)
	plt.gca().xaxis.tick_top()
	plt.gca().xaxis.set_label_position("top")
	if subplot: return
	if key in pk.SU_KEY_LIST:
		plt.xlabel(key.upper(),fontsize='large')
	else:
		plt.xlabel("Trace number",fontsize='large')
	plt.ylabel("Time (s)",fontsize='large')

def plot_wiggle(self,figsize=[5,10],fill=True,perc=100,scale=1,subplot=False):
	# plot wiggle traces
	# figsize=[5,10]: matplotlib figure size [inch]
	# fill=True: fill values greater than zero
	# perc=100: percent clip
	# scale=1: scale trace plots
	plotdata=self.data
	plotdata=perc_clip(plotdata,perc)
	print("min=%s max=%s"%(plotdata.min(),plotdata.max()))
	maxval=np.abs(plotdata).max()
	ns=pk.get_ns(self)
	dt=pk.get_dt(self)
	ntr=pk.get_ntr(self)
	t=np.arange(ns)*dt

	if not subplot: plt.figure(figsize=figsize)
	for itr in range(ntr):
		trace=plotdata[itr,:]
		x=itr+trace/maxval*scale
		plt.plot(x,t,'k-',linewidth=0.5)
		if fill: plt.fill_betweenx(t,x,itr,where=x>itr,color='black',linewidth=0.)

	plt.xlim([-2,ntr+1])
	plt.ylim([t[-1],t[0]])
	plt.gca().xaxis.tick_top()
	plt.gca().xaxis.set_label_position('top')
	if subplot: return
	plt.xlabel('Trace number',fontsize='large')
	plt.ylabel('Time (s)',fontsize='large')

def specfx(self,perc=100,cmap="jet",subplot=False):
	# Plot amplitude spectrum [DB]
	# perc=100: clip by percent
	# cmap='jet': matplotlib colormap
	data_f = np.fft.fft(perc_clip(self.data,perc),axis=1).T
	nt,nx=data_f.shape
	nf=nt/2+1
	dt=pk.get_dt(self)
	fmax=0.5/dt
	print("dt=%s, fmax=%s"%(dt,fmax))

	if not subplot: plt.figure(figsize=[5,10])
	plt.imshow(20.*np.log10(np.abs(data_f[:nf,:])),aspect="auto",extent=[0,nx-1,fmax,0],cmap=cmap)
	plt.gca().xaxis.tick_top()
	plt.gca().xaxis.set_label_position("top")
	if subplot: return
	plt.xlabel("Trace number",fontsize="large")
	plt.ylabel("Frequency (Hz)",fontsize="large")

def specfk(self,dx=0,perc=100,cmap="jet",subplot=False):
	# Plot fk spectrum [DB]
	# dx: receiver interval
	# perc=100: clip by percent
	# cmap='jet': matplotlib colormap
	data_fk = np.fft.fft2(perc_clip(self.data,perc)).T
	nt,nx=data_fk.shape
	if dx == 0:
		offset=pk.get_key(self,"offset")[:2]
		dx=np.abs(offset[1]-offset[0])/1000. # in km if offset is in meter
		if dx == 0:
			dx=1
	nf=nt/2+1
	dt=pk.get_dt(self)
	fmax=0.5/dt
	print("dt=%s, fmax=%s"%(dt,fmax))
	nk=nx/2+1
	kmax=0.5/dx
	print("dx=%s, kmax=%s"%(dx,kmax))

	if not subplot: plt.figure(figsize=[5,10])
	plt.imshow(20.*np.log10(np.fft.fftshift(np.abs(data_fk[:nf,:nk]),axes=1)),aspect="auto",extent=[-kmax,kmax,fmax,0],cmap=cmap)
	plt.gca().xaxis.tick_top()
	plt.gca().xaxis.set_label_position("top")
	if subplot: return
	plt.xlabel("Wavenumber",fontsize="large")
	plt.ylabel("Frequency (Hz)",fontsize="large")

def seis_env_dB(trc,trcgain,tnum=-1):
	D=trc.data.copy()
	Dg=trcgain.data.copy()

	ns=pk.get_ns(trc)
	dt=pk.get_dt(trc)
	t=np.arange(ns)*dt

	nx=D.shape[0]
	if tnum==-1:
		trc=np.mean(D,axis=0)
		gtrc=np.mean(Dg,axis=0)
		yl='average trace'
	elif tnum<nx:
		trc=D[tnum,:]
		gtrc=Dg[tnum,:]
		yl='trace %d'%tnum
	else:
		print 'tnum should be smaller than %d'%nx
		return

	import scipy.signal
	trc/=trc.max()
	trc_env_dB=20.*np.log10(np.abs(scipy.signal.hilbert(trc)))

	gtrc/=gtrc.max()
	gtrc_env_dB=20.*np.log10(np.abs(scipy.signal.hilbert(gtrc)))

	plt.figure()
	plt.plot(t,trc_env_dB,'b-',label='Before gain correction')
	plt.plot(t,gtrc_env_dB,'g--',label='After gain correction')
	plt.legend(loc='lower left')
	plt.grid()
	plt.xlabel('Time (s)',fontsize='large')
	plt.ylabel('Envelope (dB) of the %s'%yl,fontsize='large')

def stacking_chart(self):
	# plot source/receiver geometry
	plt.figure(figsize=[10,10])
	sx=pk.get_key(self,"sx")
	gx=pk.get_key(self,"gx")
	plt.plot(gx,sx,"ws")

	xr=gx.max()-gx.min()
	plt.xlim([gx.min()-0.1*xr,gx.max()+0.1*xr])
	yr=sx.max()-sx.min()
	plt.ylim([sx.min()-0.1*yr,sx.max()+0.2*yr])

	ntr_per_csg=pk.ntr_per_shot(self)

	for ishot in range(pk.get_nshot(self)):
		pos=sx[ishot*ntr_per_csg[ishot]]
		plt.plot(pos,pos,"r*")
		if ishot % 10 == 0:
			plt.text(gx.min()-0.05*xr,pos,"%d"%ishot,horizontalalignment="center",verticalalignment="center")

	plt.xlabel("x-axis locations",fontsize="large")
	plt.ylabel("Shot number",fontsize="large")
	plt.tick_params(axis='y',left="off",labelleft="off")
	plt.legend(["Receivers","Sources"],numpoints=1,fontsize="large")


def plot_vel(vel,h,figsize=[15,4],unit='km/s',tick=np.arange(1.5,6,1)):
	xmax=(vel.shape[0]-1)*h
	zmax=(vel.shape[1]-1)*h
	fs=14
	plt.figure(figsize=figsize)
	ax=plt.imshow(vel.transpose(),extent=(0,xmax,zmax,0))
	plt.tick_params(labelsize=fs)
	plt.xlabel('Distance (km)',fontsize=fs)
	plt.ylabel('Depth (km)',fontsize=fs)
	plt.gca().xaxis.tick_top()
	plt.gca().xaxis.set_label_position("top")
	#ax.axes.set_yticks(np.arange(0,zmax+1,1)) 

	cb=plt.colorbar(shrink=1.0,pad=0.01,aspect=10,ticks=tick)
	plt.clim([vel.min(),vel.max()])
	cb.set_label(unit,fontsize=fs)
	ct=plt.getp(cb.ax,'ymajorticklabels')
	plt.setp(ct,fontsize=fs)

def plot_mig(mig,h,figsize=[15,4]):
	xmax=(mig.shape[0]-1)*h
	zmax=(mig.shape[1]-1)*h
	fs=14
	plt.figure(figsize=figsize)
	ax=plt.imshow(mig.T,extent=(0,xmax,zmax,0),cmap=plt.cm.gray)
	plt.xlabel('Distance (km)',fontsize=fs)
	plt.ylabel('Depth (km)',fontsize=fs)
	plt.gca().xaxis.tick_top()
	plt.gca().xaxis.set_label_position("top")
	#ax.axes.set_yticks(np.arange(0,zmax+1,1))
