import cython
import numpy as np
cimport numpy as np
from . import pkbase as pk

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def impulse(int l,dtype=np.float64):
	# return impulse array with size=l
	# l: size of impulse array
	# dtype: numpy data type
	arr=np.zeros(l,dtype=dtype)
	arr[0]=1
	return arr

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def autocorr(np.ndarray[double] x):
	# auto correlation
	cdef np.ndarray[double] result = np.correlate(x, x, mode='full')
	return result[result.size/2:]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def auto_correlation_map(self,double max_lag=0.2):
	# Auto correlation of each traces
	# max_lag: maximum time lags to calculate autocorrelation [seconds]
	# output: auto-correlated SeismicTrace
	cdef double dt=pk.get_dt(self)
	cdef int N=int(np.rint(max_lag/dt))
	cdef int nx=self.data.shape[0]
	cdef int nt=self.data.shape[1]
	cdef np.ndarray[double,ndim=2] trace=np.zeros((nx,N),dtype=np.float64)
	cdef int ix
	for ix in range(nx):
		trace[ix,:]=autocorr(self.data[ix,:])[:N]
	head=self.header.copy()
	head['ns']=N
	out=pk.SeismicTrace(head,trace,self.logs())
	out.add_log("auto_correlation_map: max_lag=%s"%max_lag)
	return out

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def corr_same_len(np.ndarray[double] x,int N,np.ndarray[double] d):
	# cross correlation (x & d)
	# x, d: input arrays to correlate, same length
	# N: output array size
	cdef np.ndarray[double] rxx=np.zeros(N)
	cdef int L=len(x)
	cdef int n,l
	for n in range(N):
		for l in range(L):
			if l+n < L:
				rxx[n]+=x[l]*d[l+n]
	return rxx 

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def my_xcorr(np.ndarray[double] x,int N,np.ndarray[double] d):
	# cross correlation (x & d)
	# x, d: input arrays to correlate
	# N: output array size
	cdef int L
	if len(x)==len(d):
		return corr_same_len(x,N,d)
	else:
		L=np.max(len(d),len(x))
		# zero pad
		x=np.pad(x,(0,L-len(x)),mode='constant')
		d=np.pad(d,(0,L-len(d)),mode='constant')
		return corr_same_len(x,N,d)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def stack_corr(np.ndarray[double] strace,np.ndarray[double,ndim=2] Da,double maxlags):
	# used in surface-consistant static correction
	cdef int nx=Da.shape[0]
	cdef int nt=Da.shape[1]
	cdef np.ndarray[double,ndim=2] Dout=np.zeros_like(Da)
	cdef int ix
	cdef double cmax
	for ix in range(nx):
		ctrace=my_xcorr(strace,maxlags,Da[ix,:])
		cmax=np.argmax(ctrace)
		Dout[ix,:]=np.pad(Da[ix,cmax:nt],(0,cmax),mode='constant')
	return Dout

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def spiking_decon(self,double max_lag=0.2,double mu=0.1):
	# Spiking deconvolution
	# max_lag: maximum time lags to calculate autocorrelation [seconds]
	# mu=0.1: white noise in percent (recommended: 0.1 ~ 0.3 %)
	# output: deconvolved SeismicTrace (scaled due to the impulse on the right-hand side)
	import scipy.linalg
	cdef np.ndarray[double,ndim=2] D=self.data
	cdef double dt=pk.get_dt(self)
	cdef int N=int(np.rint(max_lag/dt))
	cdef int nx=D.shape[0]
	cdef int nt=D.shape[1]
	cdef double p_noise=mu/100.
	cdef np.ndarray[double,ndim=2] Dauto=np.zeros((nx,N),dtype=np.float64)
	cdef np.ndarray[double,ndim=2] Ds=np.zeros_like(D)
	cdef int ix
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

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def scr_static(self,int cmp_start,int cmp_end,double maxlags):
	# Surface-consistent residual static correction
	# cmp_start: first cmp number to apply the correction
	# cmp_end: last cmp number to apply the correction
	# maxlags: number of samples to use in crosscorrelation [# of samples]
	# output: static-corrected SeismicTrace
	cdef np.ndarray[double,ndim=2] Dstatic=np.zeros_like(self.data)
	cdef np.ndarray[int] cdp=pk.get_key(self,'cdp')
	cdef int cmpnum
	for cmpnum in range(cmp_start,cmp_end+1):
		#print("cmp num=",cmpnum)
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

