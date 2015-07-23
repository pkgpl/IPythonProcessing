import numpy as np
import pkbase as pk

def autocorr(x):
	# auto correlation
	result = np.correlate(x, x, mode='full')
	return result[result.size/2:]

def auto_correlation_map(self,max_lag=0.2):
	# Auto correlation of each traces
	# max_lag: maximum time lags to calculate autocorrelation [seconds]
	# output: auto-correlated SeismicTrace
	dt=pk.get_dt(self)
	N=int(np.rint(max_lag/dt))
	nx,nt=self.data.shape
	trace=np.zeros((nx,N),dtype=np.float64)
	for ix in range(nx):
		trace[ix,:]=autocorr(self.data[ix,:])[:N]
	head=self.header.copy()
	head['ns']=N
	out=pk.SeismicTrace(head,trace,self.logs())
	out.add_log("auto_correlation_map: max_lag=%s"%max_lag)
	return out


def corr_same_len(x,N,d):
	# cross correlation (x & d)
	# x, d: input arrays to correlate, same length
	# N: output array size
	rxx=np.zeros(N)
	L=len(x)
	for n in range(N):
		for l in range(L):
			if l+n < L:
				rxx[n]+=x[l]*d[l+n]
	return rxx 

def my_xcorr(x,N,d):
	# cross correlation (x & d)
	# x, d: input arrays to correlate
	# N: output array size
	if len(x)==len(d):
		return corr_same_len(x,N,d)
	else:
		L=np.max(len(d),len(x))
		# zero pad
		x=np.pad(x,(0,L-len(x)),mode='constant')
		d=np.pad(d,(0,L-len(d)),mode='constant')
		return corr_same_len(x,N,d)

def stack_corr(strace,Da,maxlags):
	# used in surface-consistant static correction
	nx,nt=Da.shape
	Dout=np.zeros_like(Da)
	for ix in range(nx):
		ctrace=my_xcorr(strace,maxlags,Da[ix,:])
		cmax=np.argmax(ctrace)
		Dout[ix,:]=np.pad(Da[ix,cmax:nt],(0,cmax),mode='constant')
	return Dout

def spiking_decon(self,max_lag=0.2,mu=0.1):
	# Spiking deconvolution
	# max_lag: maximum time lags to calculate autocorrelation [seconds]
	# mu=0.1: white noise in percent (recommended: 0.1 ~ 0.3 %)
	# output: deconvolved SeismicTrace (scaled due to the impulse on the right-hand side)
	import scipy.linalg
	D=self.data
	dt=pk.get_dt(self)
	N=int(np.rint(max_lag/dt))
	nx,nt=D.shape
	p_noise=mu/100.
	Dauto=np.zeros((nx,N),dtype=np.float64)
	Ds=np.zeros_like(D)
	print("Calculating the filter")
	for ix in range(nx):
		Dauto[ix,:]=my_xcorr(D[ix,:],N,D[ix,:]) #autocorr(D[ix,:])[:N]
	DDauto=np.sum(Dauto,0)
	DDauto[0]+=np.abs(DDauto[0])*p_noise # prewhitening
	Rxx=scipy.linalg.toeplitz(DDauto)
	Rxd=impulse(N)
	h_opt=scipy.linalg.solve(Rxx,Rxd)
	print("Applying the filter")
	for ix in range(nx):
		Ds[ix,:]=np.convolve(D[ix,:],h_opt)[:nt]
		#Ds[ix,:]=scipy.signal.fftconvolve(D[ix,:],h_opt,mode='full')[:nt] # alternative
	out=pk.SeismicTrace(self.header,Ds,self.logs())
	out.add_log("spiking_decon: max_lag=%s mu=%s"%(max_lag,mu))
	return out

def impulse(l,dtype=np.float64):
	# return impulse array with size=l
	# l: size of impulse array
	# dtype: numpy data type
	arr=np.zeros(l,dtype=dtype)
	arr[0]=1
	return arr

def scr_static(self,cmp_start,cmp_end,maxlags):
	# Surface-consistent residual static correction
	# cmp_start: first cmp number to apply the correction
	# cmp_end: last cmp number to apply the correction
	# maxlags: number of samples to use in crosscorrelation [# of samples]
	# output: static-corrected SeismicTrace
	Dstatic=np.zeros_like(self.data)
	cdp=pk.get_key(self,'cdp')
	for cmpnum in range(cmp_start,cmp_end+1):
		print "cmp num=",cmpnum
		su1=pk.window(self,'cdp',cmpnum)
		Da=su1.data
		# first stacked trace
		Da1=stack_corr(np.sum(Da,0),Da,maxlags)
		# second stacked trace
		Da2=stack_corr(np.sum(Da1,0),Da,maxlags)
		Dstatic[cdp==cmpnum,:]=Da2
	out=pk.SeismicTrace(self.header,Dstatic,self.logs())
	out.add_log('scr_static: cmp_start=%s cmp_end=%s maxlags=%s'%(cmp_start,cmp_end,maxlags))
	return out

