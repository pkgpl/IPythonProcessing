import cython
import numpy as np
cimport numpy as np

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def velan(np.ndarray[double,ndim=2] D,double dt,np.ndarray[double] h,double vmin,double vmax,int nv,int R,int L):
	cdef int it,iv,ig,ih,i1,i2
	cdef double iss
	cdef int nt=D.shape[0], nh=D.shape[1]
	cdef np.ndarray[double] v=np.linspace(vmin,vmax,nv)
	cdef np.ndarray[double] tau=np.arange(0,nt-1,R)*dt
	cdef int ntau=len(tau)
	cdef int lwin=2*L+1
	cdef np.ndarray[double,ndim=2] H=(np.hamming(lwin)*np.ones((nh,1))).T
	cdef np.ndarray[double,ndim=2] S=np.zeros((ntau,nv))
	cdef np.ndarray[double] time=np.empty(nh)
	cdef np.ndarray[double] ts=np.empty(nh)
	
	for it in range(ntau): # loop over t_0
		for iv in range(nv): # loop over vel
			for ih in range(nh):
				time[ih]=(tau[it]**2+(h[ih]/v[iv])**2)**0.5
			
			s=np.zeros((lwin,nh))
			for ig in range(-L,L): # loop over window
				for ih in range(nh):
					ts[ih]=time[ih]+ig*dt
				for ih in range(nh):
					iss=ts[ih]/dt+1.0
					i1=np.floor(iss)
					i2=i1+1
					if i1>=0 and i2<nt:
						a=iss-i1
						s[ig+L,ih]=(1.-a)*D[i1,ih]+a*D[i2,ih] # linear intpl
			s*=H
			s1=np.sum(np.sum(s,1)**2)
			s2=np.sum(np.sum(s**2))
			S[it,iv]=np.abs(s1-s2)
	return S,tau,v

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def nmo(np.ndarray[double,ndim=2] D,double dt,np.ndarray[double] h,np.ndarray[double] tnmo,np.ndarray[double] vnmo,double max_stretch):
	cdef int it,ih
	cdef int nt=D.shape[0], nh=D.shape[1]
	#nt,nh=D.shape
	cdef np.ndarray[double] mute_count=np.zeros(nt)
	# interpolate the nmo velocity
	cdef np.ndarray[double] ti=np.arange(0,nt)*dt
	cdef np.ndarray[double] vi=np.interp(ti,tnmo,vnmo)
	
	cdef np.ndarray[double,ndim=2] Dnmo=np.zeros_like(D)
	cdef np.ndarray[double] M=np.zeros(nt)
	cdef double arg, time, stretch, its, a
	cdef int it1,it2
	
	for it in range(nt):
		for ih in range(nh):
			arg=ti[it]**2+(h[ih]/vi[it])**2
			time=np.sqrt(arg)
			stretch=(time-ti[it])/(ti[it]+1.e-10)
			if stretch < max_stretch/100.:
				M[it]+=1
				its=time/dt+1
				it1=np.floor(time/dt+1)
				it2=it1+1
				a=its-it1
				if it2 < nt:
					Dnmo[it,ih]=(1.-a)*D[it1,ih]+a*D[it2,ih]
	return Dnmo,M,ti,vi

