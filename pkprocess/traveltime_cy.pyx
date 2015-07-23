# coding: utf-8
import numpy as np
import cython
cimport numpy as np

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def t_side(np.ndarray[double,ndim=2] time,np.ndarray[double,ndim=2] s,double h,int ix,int iz,double ta,double tb,double tc):
    cdef double tmp=(h*s[ix,iz])**2-0.25*(tb-tc)**2
    if tmp<0.:
        time[ix,iz]=ta+h*s[ix,iz]
    else:
        time[ix,iz]=ta+np.sqrt(tmp)
        
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def t_bottom(np.ndarray[double,ndim=2] time,np.ndarray[double,ndim=2] s,double h,int ix,int iz):
    cdef double ta=time[ix,iz-1]
    cdef double tb=time[ix-1,iz-1]
    cdef double tc=time[ix+1,iz-1]
    t_side(time,s,h,ix,iz,ta,tb,tc)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def t_top(np.ndarray[double,ndim=2] time,np.ndarray[double,ndim=2] s,double h,int ix,int iz):
    cdef double ta=time[ix,iz+1]
    cdef double tb=time[ix-1,iz+1]
    cdef double tc=time[ix+1,iz+1]
    t_side(time,s,h,ix,iz,ta,tb,tc)
    
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def t_left(np.ndarray[double,ndim=2] time,np.ndarray[double,ndim=2] s,double h,int ix,int iz):
    cdef double ta=time[ix+1,iz]
    cdef double tb=time[ix+1,iz-1]
    cdef double tc=time[ix+1,iz+1]
    t_side(time,s,h,ix,iz,ta,tb,tc)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def t_right(np.ndarray[double,ndim=2] time,np.ndarray[double,ndim=2] s,double h,int ix,int iz):
    cdef double ta=time[ix-1,iz]
    cdef double tb=time[ix-1,iz-1]
    cdef double tc=time[ix-1,iz+1]
    t_side(time,s,h,ix,iz,ta,tb,tc)



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def t_side1(np.ndarray[double,ndim=2] time,np.ndarray[double,ndim=2] s,double h,int ix,int iz,double ta):
    time[ix,iz]=ta+h*s[ix,iz]

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def t_bottom1(np.ndarray[double,ndim=2] time,np.ndarray[double,ndim=2] s,double h,int ix,int iz):
    cdef double ta=time[ix,iz-1]
    t_side1(time,s,h,ix,iz,ta)
    
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def t_top1(np.ndarray[double,ndim=2] time,np.ndarray[double,ndim=2] s,double h,int ix,int iz):
    cdef double ta=time[ix,iz+1]
    t_side1(time,s,h,ix,iz,ta)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def t_left1(np.ndarray[double,ndim=2] time,np.ndarray[double,ndim=2] s,double h,int ix,int iz):
    cdef double ta=time[ix+1,iz]
    t_side1(time,s,h,ix,iz,ta)
    
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def t_right1(np.ndarray[double,ndim=2] time,np.ndarray[double,ndim=2] s,double h,int ix,int iz):
    cdef double ta=time[ix-1,iz]
    t_side1(time,s,h,ix,iz,ta)



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def t_corner(np.ndarray[double,ndim=2] time,np.ndarray[double,ndim=2] s,double h,int ix,int iz,double ta,double tb,double tc):
    cdef double tmp=2.*(h*s[ix,iz])**2-(tb-tc)**2
    if tmp<0.:
        time[ix,iz]=np.minimum(tb,tc)+h*s[ix,iz]
    else:
        time[ix,iz]=ta+np.sqrt(tmp)
        
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def t_bottom_right(np.ndarray[double,ndim=2] time,np.ndarray[double,ndim=2] s,double h,int ix,int iz):
    cdef double ta=time[ix-1,iz-1]
    cdef double tb=time[ix-1,iz]
    cdef double tc=time[ix,iz-1]
    t_corner(time,s,h,ix,iz,ta,tb,tc)
    
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def t_bottom_left(np.ndarray[double,ndim=2] time,np.ndarray[double,ndim=2] s,double h,int ix,int iz):
    cdef double ta=time[ix+1,iz-1]
    cdef double tb=time[ix+1,iz]
    cdef double tc=time[ix,iz-1]
    t_corner(time,s,h,ix,iz,ta,tb,tc)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def t_top_right(np.ndarray[double,ndim=2] time,np.ndarray[double,ndim=2] s,double h,int ix,int iz):
    cdef double ta=time[ix-1,iz+1]
    cdef double tb=time[ix-1,iz]
    cdef double tc=time[ix,iz+1]
    t_corner(time,s,h,ix,iz,ta,tb,tc)
    
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def t_top_left(np.ndarray[double,ndim=2] time,np.ndarray[double,ndim=2] s,double h,int ix,int iz):
    cdef double ta=time[ix+1,iz+1]
    cdef double tb=time[ix+1,iz]
    cdef double tc=time[ix,iz+1]
    t_corner(time,s,h,ix,iz,ta,tb,tc)



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def traveltime(np.ndarray[double,ndim=2] vel,double h,int srcx,int srcz):
    cdef int nx=vel.shape[0],nz=vel.shape[1]
    cdef np.ndarray[double,ndim=2] s=1./vel

    cdef int tinit=2**10
    cdef np.ndarray[double,ndim=2] time=np.ones((nx,nz))*tinit

    cdef int boundary_top=srcz
    cdef int boundary_bottom=srcz
    cdef int boundary_left=srcx
    cdef int boundary_right=srcx
    cdef int iring,ix,ix_l,ix_r,iz,iz_up,iz_dn,imin

# near source
    time[srcx,srcz]=0.
    # side
    if srcz > 0: # top
        time[srcx,srcz-1]=h/2.*(s[srcx,srcz]+s[srcx,srcz-1])
        boundary_top=srcz-1
    if srcz < nz-1: # bottom
        time[srcx,srcz+1]=h/2.*(s[srcx,srcz]+s[srcx,srcz+1])
        boundary_bottom=srcz+1
    if srcx > 0: # left
        time[srcx-1,srcz]=h/2.*(s[srcx,srcz]+s[srcx-1,srcz])
        boundary_left=srcx-1
    if srcx < nx-1: # right
        time[srcx+1,srcz]=h/2.*(s[srcx,srcz]+s[srcx+1,srcz])
        boundary_right=srcx+1

    # corners
    if srcz > 0 and srcx > 0:
        t_top_left(time,s,h,srcx-1,srcz-1)
    if srcz > 0 and srcx < nx-1:
        t_top_right(time,s,h,srcx+1,srcz-1)
    if srcz < nz-1 and srcx > 0:
        t_bottom_left(time,s,h,srcx-1,srcz+1)
    if srcz < nz-1 and srcx < nx-1:
        t_bottom_right(time,s,h,srcx+1,srcz+1)

# ring        
    for iring in range(max(nx,nz)):
    # Bottom edge

        if boundary_bottom < nz-1:
            # z-index to calculate
            iz=boundary_bottom+1
            iz_up=boundary_bottom

            # find minimum traveltime
            imin=np.argmin(time[:,boundary_bottom])
            t_bottom(time,s,h,imin,iz)

            # find local min.
            for ix in range(boundary_left+1,boundary_right):
                if time[ix,iz_up]<time[ix-1,iz_up] and time[ix,iz_up]<time[ix+1,iz_up]:
                    t_bottom(time,s,h,ix,iz)

            # if marginal elements are local min.
            if time[boundary_left,iz_up]<time[boundary_left+1,iz_up]:
                t_bottom1(time,s,h,boundary_left,iz)
            if time[boundary_right,iz_up]<time[boundary_right-1,iz_up]:
                t_bottom1(time,s,h,boundary_right,iz)

            # left
            # incoming
            for ix in range(boundary_left+1,imin):
                if time[ix-1,iz_up]<time[ix,iz_up]:
                    t_bottom_right(time,s,h,ix,iz)
            # outgoing
            for ix in range(imin-1,boundary_left-1,-1):
                if time[ix+1,iz_up]<time[ix,iz_up]:
                    t_bottom_left(time,s,h,ix,iz)

            # right
            # incoming
            for ix in range(boundary_right-1,imin,-1):
                if time[ix+1,iz_up]<time[ix,iz_up]:
                    t_bottom_left(time,s,h,ix,iz)
            # outgoing
            for ix in range(imin+1,boundary_right+1):
                if time[ix-1,iz_up]<time[ix,iz_up]:
                    t_bottom_right(time,s,h,ix,iz)

            boundary_bottom+=1

    # top edge
        if boundary_top > 0:
            iz=boundary_top-1
            iz_dn=boundary_top

            imin=np.argmin(time[:,boundary_top])
            t_top(time,s,h,imin,iz)
            # local min
            for ix in range(boundary_left+1,boundary_right):
                if time[ix,iz_dn]<time[ix-1,iz_dn] and time[ix,iz_dn]<time[ix+1,iz_dn]:
                    t_top(time,s,h,ix,iz)
                    
            # if marginal elements are local min.
            if time[boundary_left,iz_dn]<time[boundary_left+1,iz_dn]:
                t_top1(time,s,h,boundary_left,iz)
            if time[boundary_right,iz_dn]<time[boundary_right-1,iz_dn]:
                t_top1(time,s,h,boundary_right,iz)

            # left
            # incoming
            for ix in range(boundary_left+1,imin):
                if time[ix-1,iz_dn]<time[ix,iz_dn]:
                    t_top_right(time,s,h,ix,iz)
            # outgoing
            for ix in range(imin-1,boundary_left-1,-1):
                if time[ix+1,iz_dn]<time[ix,iz_dn]:
                    t_top_left(time,s,h,ix,iz)

            # right
            # incoming
            for ix in range(boundary_right-1,imin,-1):
                if time[ix+1,iz_dn]<time[ix,iz_dn]:
                    t_top_left(time,s,h,ix,iz)
            # outgoing
            for ix in range(imin+1,boundary_right+1):
                if time[ix-1,iz_dn]<time[ix,iz_dn]:
                    t_top_right(time,s,h,ix,iz)

            boundary_top-=1

    # left edge
        if boundary_left > 0:
            ix=boundary_left-1
            ix_r=boundary_left

            imin=np.argmin(time[boundary_left,:])
            t_left(time,s,h,ix,imin)
            # local min
            for iz in range(boundary_top+1,boundary_bottom):
                if time[ix_r,iz]<time[ix_r,iz-1] and time[ix_r,iz]<time[ix_r,iz+1]:
                    t_left(time,s,h,ix,iz)

            # if marginal elements are local min.
            if time[ix_r,boundary_top]<time[ix_r,boundary_top+1]:
                t_left1(time,s,h,ix,boundary_top)
            if time[ix_r,boundary_bottom]<time[ix_r,boundary_bottom-1]:
                t_left1(time,s,h,ix,boundary_bottom)
                
            # top
            # incoming
            for iz in range(boundary_top+1,imin):
                if time[ix_r,iz-1]<time[ix_r,iz]:
                    t_bottom_left(time,s,h,ix,iz)
            # outgoing
            for iz in range(imin-1,boundary_top-1,-1):
                if time[ix_r,iz+1]<time[ix_r,iz]:
                    t_top_left(time,s,h,ix,iz)

            # bottom
            # incoming
            for iz in range(boundary_bottom-1,imin,-1):
                if time[ix_r,iz+1]<time[ix_r,iz]:
                    t_top_left(time,s,h,ix,iz)
            # outgoing
            for iz in range(imin+1,boundary_bottom+1):
                if time[ix_r,iz-1]<time[ix_r,iz]:
                    t_bottom_left(time,s,h,ix,iz)

            boundary_left-=1

    # right edge
        if boundary_right < nx-1:
            ix=boundary_right+1
            ix_l=boundary_right

            imin=np.argmin(time[boundary_right,:])
            t_right(time,s,h,ix,imin)
            # local min
            for iz in range(boundary_top+1,boundary_bottom):
                if time[ix_l,iz]<time[ix_l,iz-1] and time[ix_l,iz]<time[ix_l,iz+1]:
                    t_right(time,s,h,ix,iz)
                    
            # if marginal elements are local min.
            if time[ix_l,boundary_top]<time[ix_l,boundary_top+1]:
                t_right1(time,s,h,ix,boundary_top)
            if time[ix_l,boundary_bottom]<time[ix_l,boundary_bottom-1]:
                t_right1(time,s,h,ix,boundary_bottom)
                
            # top
            # incoming
            for iz in range(boundary_top+1,imin):
                if time[ix_l,iz-1]<time[ix_l,iz]:
                    t_bottom_right(time,s,h,ix,iz)
            # outgoing
            for iz in range(imin-1,boundary_top-1,-1):
                if time[ix_l,iz+1]<time[ix_l,iz]:
                    t_top_right(time,s,h,ix,iz)

            # bottom
            # incoming
            for iz in range(boundary_bottom-1,imin,-1):
                if time[ix_l,iz+1]<time[ix_l,iz]:
                    t_top_right(time,s,h,ix,iz)
            # outgoing
            for iz in range(imin+1,boundary_bottom+1):
                if time[ix_l,iz-1]<time[ix_l,iz]:
                    t_bottom_right(time,s,h,ix,iz)

            boundary_right+=1

    # corners
        ix=boundary_left
        iz=boundary_top
        if not time[ix,iz] < tinit:
            t_top_left(time,s,h,ix,iz)

        ix=boundary_right
        iz=boundary_top
        if not time[ix,iz] < tinit:
            t_top_right(time,s,h,ix,iz)

        ix=boundary_left
        iz=boundary_bottom
        if not time[ix,iz] < tinit:
            t_bottom_left(time,s,h,ix,iz)

        ix=boundary_right
        iz=boundary_bottom
        if not time[ix,iz] < tinit:
            t_bottom_right(time,s,h,ix,iz)
    
    return time


