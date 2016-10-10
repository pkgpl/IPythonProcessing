from numpy import *
from matplotlib.pyplot import *
import pickle
#import pkprocess as pk
from . import pkbase as pk

UNIT="m/s"

from __init__ import CYTHON
if CYTHON:
	from .velan_cy import velan,nmo
else:
	from .velan import velan,nmo

class Picker:
	def __init__(self,Da,dt,S,tau,v,cmpnum,h,max_stretch):
		self.fig=figure(figsize=[15,10])
		self.Da=Da
		self.dt=dt
		self.S=S
		self.tau=tau
		self.v=v
		self.cmpnum=cmpnum
		self.offset=h
		self.max_stretch=max_stretch
		self.vstack=zeros(0)
		self.tstack=zeros(0)
		self.pvstack=zeros(0)
		self.ptstack=zeros(0)

	def draw(self):
		self.ax_left=self.fig.add_subplot(1,3,1)
		self.plot_wigb_sub(self.ax_left)
		ylabel('Time (s)',fontsize='large')
		
		self.ax_vel=self.fig.add_subplot(1,3,2)
		self.plot_vel_spec(self.ax_vel)
		title("CMP: %s"%self.cmpnum,fontsize='large',y=1.08)

		self.vline,=self.ax_vel.plot(self.vstack,self.tstack,'r-*',ms=10)
		xlim([self.v[0],self.v[-1]])
		ylim([self.tau[-1],self.tau[0]])

		self.ax_right=self.fig.add_subplot(1,3,3)
		self.plot_wigb_sub(self.ax_right)
		show()
		print('num picks=',len(self.vstack))
		return self.vstack,self.tstack

	def redraw(self):
		# vel pick line - center panel
		self.vline.set_xdata(self.vstack)
		self.vline.set_ydata(self.tstack)
		# nmo cmp - right panel
		if len(self.vstack)>0:
			for itr in range(self.ntr):
				trace=self.Dnmo[:,itr]
				self.nmo_cmp_line[itr].set_xdata(itr+trace/self.maxval)
		self.fig.canvas.draw()

	def plot_wigb_sub(self,ax):
		D=self.Da
		dt=self.dt
		maxval=D.max()/2.

		nt,ntr=D.shape
		self.ntr=ntr
		self.maxval=maxval
		self.nmo_cmp_line=range(ntr)
		t=arange(nt)*dt
		for itr in range(ntr):
			trace=D[:,itr]
			x=itr+trace/maxval
			line,=ax.plot(x,t,'k-',linewidth=0.5)
			#ax.fill_betweenx(t,x,itr,where=x>itr,color='black',linewidth=0.)
			self.nmo_cmp_line[itr]=line

		xlim([-2,ntr+1])
		ylim([t[-1],t[0]])
		xlabel('Trace number',fontsize='large')
		ax.xaxis.tick_top()
		ax.xaxis.set_label_position('top')

	def plot_vel_spec(self,ax):
		S=self.S
		tau=self.tau
		v=self.v
		im=ax.imshow(S,extent=[v[0],v[-1],tau[-1],tau[0]],aspect='auto',cmap='jet')
		grid()
		colorbar(im,ax=ax)
		xlabel('Velocity (%s)'%UNIT,fontsize='large')
		ax.xaxis.tick_top()
		ax.xaxis.set_label_position('top')

	def onclick(self,event):
		if event.inaxes != self.ax_vel:
			return
		if event.button==2:
			return
		if event.button==1: # left click to pick
			self.vstack=append(self.vstack,event.xdata)
			self.tstack=append(self.tstack,event.ydata)
			self.pvstack=append(self.pvstack,event.x)
			self.ptstack=append(self.ptstack,event.y)
			# sort by time
			indx=argsort(self.tstack)
			self.tstack=self.tstack[indx]
			self.vstack=self.vstack[indx]
			self.ptstack=self.ptstack[indx]
			self.pvstack=self.pvstack[indx]
			print('xdata=%s, ydata=%s, npick=%s'%(event.xdata, event.ydata, len(self.vstack)))
			#print('x=%s, y=%s'%(event.x,event.y))
			if len(self.vstack)>0:
				self.Dnmo,M,ti,vi=nmo(self.Da,self.dt,self.offset.astype(np.float64),self.tstack,self.vstack,float(self.max_stretch))
		if event.button==3: # right click to remove
			if len(self.vstack)>0:
				self.remove(event.x,event.y)
			if len(self.vstack)>0:
				self.Dnmo,M,ti,vi=nmo(self.Da,self.dt,self.offset.astype(np.float64),self.tstack,self.vstack,float(self.max_stretch))
		self.redraw()
	
	def remove(self,x,y):
		# remove nearest y(=time)
		i=argmin(sqrt((self.pvstack-x)**2+(self.ptstack-y)**2))
		self.tstack=delete(self.tstack,i)
		self.vstack=delete(self.vstack,i)
		self.ptstack=delete(self.ptstack,i)
		self.pvstack=delete(self.pvstack,i)


def vt_interp(vtlist,cmpn,cmp_num):
	if cmpn in cmp_num:
		return vtlist[np.argwhere(cmp_num==cmpn)]
	else:
		cmp_left=cmp_num[cmp_num<cmpn].max()
		cmp_right=cmp_num[cmp_num>cmpn].min()
		wl=float(cmp_right-cmpn)/float(cmp_right-cmp_left)
		wr=1.-wl

		vleft,tleft=vtlist[np.argwhere(cmp_num==cmp_left)]
		vright,tright=vtlist[np.argwhere(cmp_num==cmp_right)]
		tnmo=np.unique(np.hstack((tleft,tright)))
		vleft_interp=np.interp(tnmo,tleft,vleft)
		vright_interp=np.interp(tnmo,tright,vright)
		vnmo=vleft_interp*wl+vright_interp*wr
		return vnmo,tnmo

def vel_picking_nmo(su,vmin,dv,nv,cmp_start,cmp_end,cmp_step,max_stretch,R,L):
	Dsort=su.data.T
	cmps=pk.get_key(su,'cdp')
	dt=pk.get_dt(su)
	#ns=pk.get_ns(su)

	#t=np.arange(ns)*dt
	vmax=(nv-1)*dv+vmin
	cmp_num=np.arange(cmp_start,cmp_end+cmp_step,cmp_step)
	vtlist=range(len(cmp_num))

	for icmp,cmpnum in enumerate(cmp_num):
		su1=pk.window(su,'cdp',cmpnum)
		Da=su1.data.T
		h=pk.get_key(su1,'offset')
		nt,nx=Da.shape
		print("CMP number: %s, num folds: %s"%(cmpnum,nx))
		print("Calculating the semblance panel...")
		S,tau,v=velan(Da,dt,h.astype(np.float64),vmin,vmax,nv,R,L)
		if nx<2:
			print("pass:: num folds = %s"%nx)
			continue
		
		picker=Picker(Da,dt,S,tau,v,cmpnum,h,max_stretch)
		cid=picker.fig.canvas.mpl_connect('button_press_event',picker.onclick)
		vnmo,tnmo=picker.draw()
		vtlist[icmp]=(vnmo,tnmo)
	with open("picks.pickle","wb") as f:
		pickle.dump(vtlist,f,pickle.HIGHEST_PROTOCOL)
	#with open("picks.pickle","rb") as f:
	#	vtlist=pickle.load(f)
	
	Dout=np.zeros_like(Dsort)
	for icmp,cmpnum in enumerate(cmp_num):
		# Normal MoveOut
		if cmpnum < cmp_end:
			for ll in range(cmp_step):
				cmpn=cmpnum+ll
				print("applying nmo to cmp #%s"%cmpn)

				# interpolate nmo velocity
				vnmo,tnmo=vt_interp(vtlist,cmpn,cmp_num)

				su1=pk.window(su,'cdp',cmpn)
				Da=su1.data.T
				h=pk.get_key(su1,'offset')
				Dnmo,M,ti,vi=nmo(Da,dt,h.astype(np.float64),tnmo,vnmo,float(max_stretch))
				jj=(cmps == cmpn)
				Dout[:,jj]=Dnmo
	out=pk.SeismicTrace(su.header,Dout.T,su.logs())
	#out.nmo_picks=vtlist
	out.nmo_picks={cmpnum:vt for (cmpnum, vt) in zip(cmp_num,vtlist)} # dictionary
	out.add_log('velan_nmo: vmin=%s dv=%s nv=%s cmp_start=%s cmp_end=%s cmp_step=%s max_stretch=%s R=%s L=%s'%(vmin,dv,nv,cmp_start,cmp_end,cmp_step,max_stretch,R,L))
	return out
		

#input_filename="marm_gain_decon_agc_sort.trc"
#output_filename="marm_gain_decon_agc_sort_nmo.trc"
#
#vmin=1500.
#dv=100
#nv=41
#
#cmp_start=25
#cmp_end=550
#cmp_step=25
#
#max_stretch=10 # max. allowable nmo stretch in %
#
#R = 4 # 2
#L = 8 # 10
#
#su = pk.from_trace(input_filename)
#
#so = vel_picking_nmo(su,vmin,dv,nv,cmp_start,cmp_end,cmp_step,max_stretch,R,L)
#
#print("saving file: %s"%output_filename)
#so.dump(output_filename)
