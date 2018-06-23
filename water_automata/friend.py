import numpy as np

def warp(z):
    w = np.zeros((z.shape[0]+2,z.shape[1]+2,*z.shape[2:]))
    w[1:-1,1:-1] = z.copy()
    w[0,:] = w[-2,:]
    w[-1,:] = w[1,:]
    w[:,0] = w[:,-2]
    w[:,-1] = w[:,1]
    return w

def spread4(x,y):
    # Creates a matrix of 4 positive components for each cell.
    assert(x.shape==y.shape)
    assert(len(x.shape)==2)
    ds = np.zeros((*x.shape,4))
    ds[:,:,0] = np.maximum(0,x)
    ds[:,:,1] = -np.minimum(0,x)
    ds[:,:,2] = np.maximum(0,y)
    ds[:,:,3] = -np.minimum(0,y)
    return ds

def central_difference(z):
    z = warp(z)
    cd_x = z[1:-1,2:]-z[1:-1,:-2]
    cd_y = z[2:,1:-1]-z[:-2,1:-1]
    return cd_x,cd_y


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

def transference(trnf):
    # Compute the change for each individual cell, from a matrix that has
    # 4 positive components, one for each direction, per cell.
    assert(len(trnf.shape)==3)
    delta = np.zeros(trnf.shape[:2])
    delta -= np.sum(trnf,axis=2)
    delta += warp_shift_to_dir(trnf[:,:,0],0)
    delta += warp_shift_to_dir(trnf[:,:,1],1)
    delta += warp_shift_to_dir(trnf[:,:,2],2)
    delta += warp_shift_to_dir(trnf[:,:,3],3)
    return delta

def surrounding_advantage(z):
    w = warp(z)
    min_x = np.minimum(w[1:-1,2:],w[1:-1,:-2])
    min_y = np.minimum(w[2:,1:-1],w[:-2,1:-1])
    adv_x = np.maximum(0,z - min_x)
    adv_y = np.maximum(0,z - min_y)
    return adv_x,adv_y

    # dx = np.maximum(z-np.maximum(w[1:-1,2:],w[1:-1,:-2]),0)
    # dy = np.maximum(z-np.maximum(w[2:,1:-1],w[:-2,1:-1]),0)
    # ds = np.zeros((*z.shape,4))
    # ds[:,:,0] = dx/2.0
    # ds[:,:,1] = dx/2.0
    # ds[:,:,2] = dy/2.0
    # ds[:,:,3] = dy/2.0
    return ds


# NOTE: Momentum goes divided by the density.

def friend_transfer_step(wa,mo_x,mo_y,te,delta_t=0.01,hh=1.0,gg=9.8,mo_loss_due_te=0.0):
    assert(mo_x.shape==wa.shape)
    assert(mo_y.shape==wa.shape)
    # Level and auxs
    le = te + wa
    te_l = warp(te)
    le_l = warp(le)
    #
    wa4 = np.stack([wa]*4,axis=2)
    wa4_div = wa4.copy()
    wa4_div[wa4_div==0] = np.inf
    # Contact with water or air
    cw_d = np.maximum(.0,np.minimum(wa,le_l[1:-1,1:-1]-te_l[1:-1,2:]))
    cw_a = np.maximum(.0,np.minimum(wa,le_l[1:-1,1:-1]-te_l[1:-1,:-2]))
    cw_s = np.maximum(.0,np.minimum(wa,le_l[1:-1,1:-1]-te_l[2:,1:-1]))
    cw_w = np.maximum(.0,np.minimum(wa,le_l[1:-1,1:-1]-te_l[:-2,1:-1]))
    cw = np.stack([cw_d,cw_a,cw_s,cw_w],axis=2)
    # Contact with terrain: NOTE: may be useful.
    ct = wa4 - cw
    # Contact coefficient:
    cont_coef = np.minimum(1.0,cw/wa4_div)
    # Discharges:
    mo4 = spread4(mo_x,mo_y)
    ideal_disch = mo4*delta_t/hh # ideal discharge (without terrain)
    disch = cont_coef*ideal_disch
    loss_disch = ideal_disch - disch # loss discharge (due terrain)
    n_wa = wa + transference(disch)
    # momentum transfer due discharge:
    mo_x4 = np.stack([mo_x]*4,axis=2)
    mo_y4 = np.stack([mo_y]*4,axis=2)
    mo_x_t = mo_x4*disch/wa4_div
    mo_y_t = mo_y4*disch/wa4_div
    # momentum generation due potential energy change
    cd_x,cd_y = central_difference(le)
    cd_x /= (2.0*hh)
    cd_y /= (2.0*hh)
    adv_x,adv_y = surrounding_advantage(wa)
    gen_mo_x = -0.5*cd_x*gg*delta_t*adv_x
    gen_mo_y = -0.5*cd_y*gg*delta_t*adv_y
    # momentum loss due terrain:
    loss_mo_x = np.sum(mo_x4*loss_disch/wa4_div,axis=2)*mo_loss_due_te
    loss_mo_y = np.sum(mo_y4*loss_disch/wa4_div,axis=2)*mo_loss_due_te
    # update momentums
    n_mo_x = mo_x + transference(mo_x_t) + gen_mo_x - loss_mo_x
    n_mo_y = mo_y + transference(mo_y_t) + gen_mo_y - loss_mo_y
    # Return:
    return n_wa,n_mo_x,n_mo_y,(n_mo_x,n_mo_y)

if __name__ == "__main__":
    siz = (5,3)
    #
    te = np.zeros(siz)
    mo_x = np.zeros(siz)
    mo_y = np.zeros(siz)
    # water:
    wa = np.zeros(siz)
    wa[1:3,1:3] = 5
    #
    for t in range(1000):
        print("--- t=%d:"%t)
        print("wa:")
        print(wa)
        print(np.sum(wa))
        print("mo_x:")
        print(mo_x)
        print("mo_y:")
        print(mo_y)
        # Update
        wa,mo_x,mo_y,_ = friend_transfer_step(wa,mo_x,mo_y,te)
        input()
