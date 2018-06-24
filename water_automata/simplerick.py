import numpy as np

# Sixteen iterations off the inundation model, there’s a model that works more with discrete updates than diferential equations.

# His name is Simple Rick, but it's not dummy.

# He realized long ago that the greatest thing he’d ever create were his water discharges.

# We captured that moment and run it on a loop through your RAM.

# And the floats that makes this simulation work go into every simple output matrix.

# Come home to the impossible flavor of your cc4 completion. Come home to Simple Rick’s.

def warp(z):
    w = np.zeros((z.shape[0]+2,z.shape[1]+2,*z.shape[2:]))
    w[1:-1,1:-1] = z
    w[0,:] = w[-2,:]
    w[-1,:] = w[1,:]
    w[:,0] = w[:,-2]
    w[:,-1] = w[:,1]
    return w

def warp_shift_to_dir(z,d):
    # Shift the matrix on the given direction:
    w = np.zeros(z.shape)
    if d==0:
        w[:,1:] = z[:,:-1]
        w[:,0] = z[:,-1]
    elif d==1:
        w[:,:-1] = z[:,1:]
        w[:,-1] = z[:,0]
    elif d==2:
        w[1:,:] = z[:-1,:]
        w[0,:] = z[-1,:]
    elif d==3:
        w[:-1,:] = z[1:,:]
        w[-1,:] = z[0,:]
    return w

def derivate(z):
    over = np.zeros((z.shape[0],z.shape[1],4))
    zz = warp(z)
    over[:,:,0] = zz[1:-1,2:] - z
    over[:,:,1] = z - zz[1:-1,:-2]
    over[:,:,2] = zz[2:,1:-1] - z
    over[:,:,3] = z - zz[:-2,1:-1]
    return over

def difference(z):
    over = np.zeros((z.shape[0],z.shape[1],4))
    zz = warp(z)
    over[:,:,0] = z - zz[1:-1,2:]
    over[:,:,1] = z - zz[1:-1,:-2]
    over[:,:,2] = z - zz[2:,1:-1]
    over[:,:,3] = z - zz[:-2,1:-1]
    return over

def transference(trnf,discount_original):
    # Compute the change for each individual cell, from a matrix that has
    # 4 positive components, one for each direction, per cell.
    assert(len(trnf.shape)==3)
    delta = np.zeros(trnf.shape[:2])
    if discount_original:
        delta -= np.sum(trnf,axis=2)
    delta += warp_shift_to_dir(trnf[:,:,0],0)
    delta += warp_shift_to_dir(trnf[:,:,1],1)
    delta += warp_shift_to_dir(trnf[:,:,2],2)
    delta += warp_shift_to_dir(trnf[:,:,3],3)
    return delta

def central_difference(z):
    z = warp(z)
    cd_x = z[1:-1,2:]-z[1:-1,:-2]
    cd_y = z[2:,1:-1]-z[:-2,1:-1]
    return cd_x/2.0,cd_y/2.0

def simplerick_step(wa,mo_x,mo_y,delta_t=0.01,hh=1.0,gg=9.8,vel=1.0):
    assert(mo_x.shape==wa.shape)
    assert(mo_y.shape==wa.shape)
    # Level and auxs
    le = wa
    stack = lambda v : np.stack([v]*4,axis=2)
    wa4_div = stack(wa)
    wa4_div[wa4_div==0] = np.inf
    # Native discharges:
    disch = np.maximum(difference(le),0)*vel*delta_t/hh
    # Discharges due momentum:
    disch[:,:,0] += np.maximum(0,mo_x) * delta_t/hh
    disch[:,:,1] += np.maximum(0,-mo_x) * delta_t/hh
    disch[:,:,2] += np.maximum(0,mo_y) * delta_t/hh
    disch[:,:,3] += np.maximum(0,-mo_y) * delta_t/hh
    # New water:
    n_wa = wa.copy()
    n_wa += transference(disch,True)
    # New momentums:
    n_mo_x = mo_x.copy()
    n_mo_y = mo_y.copy()
    # Momentum generated on discharges
    cd_x,cd_y = central_difference(le)
    mo_x_g = -disch*stack((2*gg*np.abs(cd_x))**0.5*np.sign(cd_x))
    mo_y_g = -disch*stack((2*gg*np.abs(cd_y))**0.5*np.sign(cd_y))
    # Momentum transfer on discharges:
    mo_x_t = stack(mo_x)*disch/wa4_div
    mo_y_t = stack(mo_y)*disch/wa4_div
    #
    n_mo_x += transference(mo_x_t,True)+transference(mo_x_g,False)
    n_mo_y += transference(mo_y_t,True)+transference(mo_y_g,False)
    # Return
    return n_wa,n_mo_x,n_mo_y

# def advantage(z):
#     dz = difference(z)
#     adv_x = np.sum(np.maximum(dz[:,:,:2],0),axis=2)
#     adv_y = np.sum(np.maximum(dz[:,:,2:],0),axis=2)
#     return adv_x,adv_y

# # Momentum generated on column
# adv_x,adv_y = advantage(le)
# mo_x_k = -0.5*gg*adv_x*cd_x*delta_t/hh
# mo_y_k = -0.5*gg*adv_y*cd_y*delta_t/hh
