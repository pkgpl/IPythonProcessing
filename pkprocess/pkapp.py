import numpy as np
from numba import jit
import scipy.signal
import scipy.interpolate
from .pkbase import *

@jit
def triang(L):
	# generate triangle
	w=np.zeros(L)
	if L%2==0: # even L
		for i in range(int(L/2)):
			n=i+1
			w[i]=(2.*n-1.)/L
		for i in range(int(L/2),L):
			n=i+1
			w[i]=2.-(2.*n-1.)/L
	else: # odd L
		for i in range(int((L+1)/2)):
			n=i+1
			w[i]=2.*n/(L+1.)
		for i in range(int((L+1)/2),L):
			n=i+1
			w[i]=2.-2.*n/(L+1.)
	return w

@jit
def gain(self,tpow=0,epow=0,agc=False,agc_gate=0.5,norm="rms"):
	# Apply gain function
	# tpow=0: t**tpow
	# epow=0: exp(epow*t)
	# agc=False: use automatic gain control (ignore tpow & epow)
	# agc_gate: agc window size [seconds]
	# norm='rms': normalize agc result by 'rms' or 'amplitude'
	# output: gained SeismicTrace
	trace = np.zeros_like(self.data)
	data = self.data
	ntr=get_ntr(self)
	ns=get_ns(self)
	dt=get_dt(self)
	if not agc:
		t = np.arange(ns)*dt
		t = t**tpow * np.exp(t*epow)
		for itr in range(ntr):
			trace[itr,:] = data[itr,:]*t
	else: # agc gain
		L=agc_gate/dt+1
		L=int(np.floor(L/2))
		h=triang(2*L+1)
		for k in range(ntr):
			e=data[k,:]**2
			rms=np.sqrt(np.convolve(e,h,'same'))
			epsi=1.e-10*np.max(rms)
			if epsi==0.: continue
			op=rms/(rms**2+epsi)
			trace[k,:]=data[k,:]*op
			if norm=='amplitude': # normalize by amplitude
				trace[k,:]/=np.max(np.abs(trace[k,:]))
			elif norm=='rms':
				trace[k,:]/=np.sqrt(np.sum(trace[k,:]**2)/ns)
	out=SeismicTrace(self.header,trace,self.logs(),self.nmo_picks)
	out.add_log("gain: tpow=%s epow=%s agc=%s agc_gate=%s norm=%s"%(tpow,epow,agc,agc_gate,norm))
	return out

@jit
def bpfilter(self,cut_off):
	# Band-pass filter
	# cut_off: [min.freq, max.freq]: frequency range to pass
	# output: band-pass filtered SeismicTrace
	dt=get_dt(self)
	nyq=0.5/dt
	b,a=scipy.signal.butter(5,np.array(cut_off)/nyq,btype='band')
	#w,h=scipy.signal.freqz(b,a)
	trace=scipy.signal.lfilter(b,a,self.data,axis=1)
	out=SeismicTrace(self.header,trace,self.logs(),self.nmo_picks)
	out.add_log("bpfilter: %s"%cut_off)
	return out


def stack(self):
	# Stack NMO-corrected CMP gathers
	# output: stacked SeismicTrace
	cmps=get_key(self,'cdp')
	cmpu=get_key_unique(self,'cdp')
	ns=get_ns(self)
	dt=get_dt(self)
	ncdp=len(cmpu)
	stacked=np.zeros((ncdp,ns))
	#for i,icmp in enumerate(cmpu):
	#	su1=window(self,'cdp',icmp)
	#	stacked[i,:]=np.sum(su1.data,axis=0)
	for i,su1 in enumerate(trace_split(self,'cdp')):
		stacked[i,:]=np.sum(su1.data,axis=0)
	head=np.zeros(stacked.shape[0],dtype=SU_HEADER_DTYPE)
	head['ns']=np.ones(stacked.shape[0],dtype=np.int32)*ns
	head['dt']=np.ones(stacked.shape[0],dtype=np.int32)*dt*1000000
	head['cdp']=cmpu
	fold_num = np.array([sum(icmp==cmps) for icmp in cmpu])
	head['shortpad']=fold_num
	out=SeismicTrace(head,stacked,self.logs(),self.nmo_picks)
	out.add_log('stack')
	return out


@jit
def stolt_mig(self,v,dx):
	# Stolt migration of CMP stacked data
	# v: constant velocity
	# dx: CMP interval
	# output: migrated SeismicTrace

	# python port of ezfkmig from http://www.biomecardio.com
	Dstacked=self.data.T
	nt,ncdp=Dstacked.shape

	dt=get_dt(self)
	num_f_pts=nt
	num_pts=num_f_pts

	U_w_kx=np.fft.fftshift(np.fft.fft2(Dstacked,(num_f_pts,num_pts)))
	# linear interpolation
	omega=2.*np.pi*np.linspace(-0.5,0.5,num_f_pts)/dt
	kx=2.*np.pi*np.linspace(-0.5,0.5,num_pts)/dx
	vv=v/np.sqrt(2.)
	kz=vv*np.sign(omega)*np.sqrt(kx**2+omega**2/vv**2)

	func=scipy.interpolate.interp2d(omega,kx,np.real(U_w_kx))
	ifunc=scipy.interpolate.interp2d(omega,kx,np.imag(U_w_kx))
	U_kz_kx=func(kz,kx)+ifunc(kz,kx)*1.0j

	Dmigrated=np.real(np.fft.ifft2(np.fft.ifftshift(U_kz_kx)))[:,:ncdp]
	out=SeismicTrace(self.header,Dmigrated.T,self.logs(),self.nmo_picks)
	out.add_log('stold_mig: v=%s dx=%s'%(v,dx))
	return out


@jit
def kirchhoff1(image,gather,times,isx,igx,dt,tdelay):
	nx=image.shape[0]
	nz=image.shape[1]
	ntr=gather.shape[0]
	nt=gather.shape[1]
	#cdef int ix,iz,itr,it
	#cdef double ts,tg,amp
	for itr in range(ntr):
		for ix in range(nx):
			for iz in range(nz):
				ts=times[isx,ix,iz]
				tg=times[igx[itr],ix,iz]
				it=int((ts+tg+tdelay)/dt)
				if it<nt:
					amp=gather[itr,it]
					image[ix,iz]+=amp
	return image


@jit
def kirchhoff(sd,h,times,tdelay):
	nx,nz=times[0].shape
	image=np.zeros((nx,nz))
	nt=get_ns(sd)
	dt=get_dt(sd)
	h_in_meter=h*1000.

	gathers=trace_split(sd,"sx")
	nshot=len(gathers)

	for ishot,gather in enumerate(gathers):
		if ishot %10 ==0:
			print(ishot,nshot)
		sx=get_key(gather,"sx")[0]
		gx=np.array(get_key(gather,"gx"))
		isx=int(sx/h_in_meter)
		igx=(gx/h_in_meter).astype(np.int32)
		image=kirchhoff1(image,gather.data,times,isx,igx,dt,tdelay)
	return image

@jit
def moving_average2d(vel,r1,r2):
	n1,n2=vel.shape
	svel=np.empty_like(vel)
	for i in range(n1):
		for j in range(n2):
			svel[i,j]=np.average(vel[max(0,i-r1):min(i+r1,n1),max(0,j-r2):min(j+r2,n2)])
	return svel

@jit
def zdiff2(img):
	dimg=np.zeros_like(img)
	nz=img.shape[1]
	for iz in range(1,nz-1):
		dimg[:,iz]=img[:,iz-1]-2.*img[:,iz]+img[:,iz+1]
	return dimg

@jit
def rmsvel(sd):
	dt=get_dt(sd)
	ns=get_ns(sd)
	at=np.array(range(ns))*dt
	dic=sd.nmo_picks
	if len(dic)==0:
		print("Please run this after velocity analysis!!")
		return
	ncmp=len(dic.keys())
	v1=np.empty((ns,ncmp))

	cmpnums=np.sort(dic.keys())
	for icmp,cmpnum in enumerate(cmpnums):
		vt=dic[cmpnum]
		v=vt[0]
		t=vt[1]
		vinterp=np.interp(at,t,v)
		v1[:,icmp]=vinterp

	cmpmin=cmpnums.min()
	cmpmax=cmpnums.max()
	cmps=get_key_unique(sd,'cdp')
	cmprange=[cmpn for cmpn in cmps if cmpmin<=cmpn and cmpn<=cmpmax]

	vrms=np.empty((ns,len(cmprange)))
	for it in range(ns):
		vrms[it,:]=np.interp(cmprange,cmpnums,v1[it,:])
	return vrms

@jit
def intervalvel(sd):
	vrms=rmsvel(sd)
	vmin=vrms.min()
	vmax=vrms.max()
	vrms=moving_average2d(vrms,50,20)
	dt=get_dt(sd)
	ns=get_ns(sd)
	at=np.array(range(ns))*dt
	vint=np.empty_like(vrms)
	vint[0,:]=vrms[0,:]
	for it in range(1,ns):
		vint[it,:]=sqrt((vrms[it,:]**2*at[it]-vrms[it-1,:]**2*at[it-1])/(at[it]-at[it-1]))
	return np.clip(vint,vmin,vmax)
