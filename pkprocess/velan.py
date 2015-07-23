import numpy as np

def velan(D,dt,h,vmin,vmax,nv,R,L):
	nt,nh=D.shape
	v=np.linspace(vmin,vmax,nv)
	
	tau=np.arange(0,nt-1,R)*dt
	ntau=len(tau)
	
	lwin=2*L+1
	
	taper=np.hamming(lwin)
	H=(taper*np.ones((nh,1))).T
	
	S=np.zeros((ntau,nv))
	for it in range(ntau): # loop over t_0
		for iv in range(nv): # loop over vel
			time=np.sqrt(tau[it]**2+(h/v[iv])**2)
			
			s=np.zeros((lwin,nh))
			for ig in range(-L,L):
				ts=time+ig*dt
				for ih in range(nh):
					iss=ts[ih]/dt+1
					i1=np.floor(iss)
					i2=i1+1
					if i1>0 and i2<nt:
						a=iss-i1
						s[ig+L,ih]=(1.-a)*D[i1,ih]+a*D[i2,ih] # linear intpl
			s*=H
			s1=np.sum(np.sum(s,1)**2)
			s2=np.sum(np.sum(s**2))
			S[it,iv]=np.abs(s1-s2)
	return S,tau,v


def nmo(D,dt,h,tnmo,vnmo,max_stretch):
	nt,nh=D.shape
	mute_count=np.zeros(nt)
	# interpolate the nmo velocity
	ti=np.arange(0,nt)*dt
	vi=np.interp(ti,tnmo,vnmo)

	Dnmo=np.zeros_like(D)
	M=np.zeros(nt)
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

