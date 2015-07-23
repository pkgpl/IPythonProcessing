import numpy as np
import cython
cimport numpy as np

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def kirchhoff1(np.ndarray[double,ndim=2] image,np.ndarray[double,ndim=2] gather,np.ndarray[float,ndim=3] times, int isx,np.ndarray[int] igx,double dt,double tdelay):
	cdef int nx=image.shape[0],nz=image.shape[1],ntr=gather.shape[0],nt=gather.shape[1]
	cdef int ix,iz,itr,it
	cdef double ts,tg,amp
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

