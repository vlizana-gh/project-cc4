import numpy as np

def warp(z):
    w = np.zeros((z.shape[0]+2,z.shape[1]+2,*z.shape[2:]))
    w[1:-1,1:-1] = z.copy()
    w[0,:] = w[-2,:]
    w[-1,:] = w[1,:]
    w[:,0] = w[:,-2]
    w[:,-1] = w[:,1]
    return w

def warp_delta(z1,z2):
    assert(z1.shape==z2.shape)
    assert(len(z1.shape)==2)
    d = np.zeros((*z1.shape,4))
    # D
    d[:,:-1,0] = z1[:,:-1] - z2[:,1:]
    d[:,-1,0] = z1[:,-1] - z2[:,0]
    # A
    d[:,1:,1] = z1[:,1:] - z2[:,:-1]
    d[:,0,1] = z1[:,0] - z2[:,-1]
    # S
    d[:-1,:,2] = z1[:-1,:] - z2[1:,:]
    d[-1,:,2] = z1[-1,:] - z2[0,:]
    # W
    d[1:,:,3] = z1[1:,:] - z2[:-1,:]
    d[0,:,3] = z1[0,:] - z2[-1,:]
    #
    return d


def warp_focus_sum(z):
    w = np.zeros(z.shape[:2])
    # D
    w[:,1:] += z[:,:-1,0]
    w[:,0] += z[:,-1,0]
    # A
    w[:,:-1] += z[:,1:,1]
    w[:,-1] += z[:,0,1]
    # S
    w[1:,:] += z[:-1,:,2]
    w[0,:] += z[-1,:,2]
    # W
    w[:-1,:] += z[1:,:,3]
    w[-1,:] += z[0,:,3]
    #
    return w

def warp_transf_on_dir(z):
    wz = warp(z)
    t = np.zeros(z.shape)
    t[:,:,0] = wz[1:-1,:-2,0]
    t[:,:,1] = wz[1:-1,2:,1]
    t[:,:,2] = wz[:-2,1:-1,2]
    t[:,:,3] = wz[2:,1:-1,3]
    return t

def bounce(z):
    assert(z.shape[2]==4)
    bz = np.zeros(z.shape)
    bz[:,:,0] = z[:,:,1]
    bz[:,:,1] = z[:,:,0]
    bz[:,:,2] = z[:,:,3]
    bz[:,:,3] = z[:,:,2]
    return bz


def energy_transfer_step(wa,kine,te,gg=9.8,delta_t=1.0,hh=1.0,nn=0.01,pp=997,watransf=0.2):
    assert(kine.shape[:-1]==wa.shape)
    assert(te.shape==wa.shape)
    assert(kine.shape[2]==4) # WASD
    # Level and auxs
    wa_l = warp(wa)
    te_l = warp(te)
    le_l = te_l+wa_l
    wa4 = np.stack([wa,wa,wa,wa],axis=2)
    wa4_div = wa4.copy()
    wa4_div[wa4_div==0] = 1
    # Contact with air
    ca_d = np.maximum(0,np.minimum(wa,le_l[1:-1,1:-1]-le_l[1:-1,2:]))
    ca_a = np.maximum(0,np.minimum(wa,le_l[1:-1,1:-1]-le_l[1:-1,:-2]))
    ca_s = np.maximum(0,np.minimum(wa,le_l[1:-1,1:-1]-le_l[2:,1:-1]))
    ca_w = np.maximum(0,np.minimum(wa,le_l[1:-1,1:-1]-le_l[:-2,1:-1]))
    ca = np.stack([ca_d,ca_a,ca_s,ca_w],axis=2)
    # Contact with water
    cw_d = np.maximum(0,np.minimum(wa,le_l[1:-1,1:-1]-te_l[1:-1,2:]))
    cw_a = np.maximum(0,np.minimum(wa,le_l[1:-1,1:-1]-te_l[1:-1,:-2]))
    cw_s = np.maximum(0,np.minimum(wa,le_l[1:-1,1:-1]-te_l[2:,1:-1]))
    cw_w = np.maximum(0,np.minimum(wa,le_l[1:-1,1:-1]-te_l[:-2,1:-1]))
    cw = np.stack([cw_d,cw_a,cw_s,cw_w],axis=2) - ca
    # Contact with terrain
    ct = wa4 - cw - ca
    # Contact with terrain coefficient:
    ct_coef = np.minimum(1,ct/wa4_div)
    # Discharge velocity
    divel = np.minimum((1.0/nn)*wa4**(2/3.0)*(ca/hh)**0.5,(wa4*gg)**0.5)
    # Height transfer due discharge
    disch_w_t = delta_t*ca*divel/hh
    # Movement velocity
    mass_div = hh*hh*pp*wa4
    mass_div[kine==0] = 1.0
    kivel = (2*kine/mass_div)**0.5
    # Height transfer due kinetic
    kine_w_t = (ca+watransf*cw)*kivel*delta_t/hh
    # Energy bounce:
    kine_e_b = (ct_coef)*kine*kivel*delta_t/hh
    # Energy transfer:
    kine_e_t = (1.0-ct_coef)*kine*kivel*delta_t/hh
    # @@@ HEIGHT TRANSFER @@@
    wa_t = np.minimum(disch_w_t+kine_w_t,wa4/5)
    n_wa = wa.copy()
    n_wa += warp_focus_sum(wa_t) - np.sum(wa_t,axis=2)
    # Kinetic energy change due change in potential energy:
    kine_new = gg*hh*hh*pp*warp_delta(wa,n_wa)*wa_t
    # New kinetic energy:
    n_kine = kine + warp_transf_on_dir(kine_new+kine_e_t+bounce(kine_e_b))
    n_kine -= kine_e_t+kine_e_b
    n_kine = np.maximum(0,n_kine)
    #
    return n_wa,n_kine

# nn = Manning coefficient of roughness - ranging from 0.01 (a clean and smooth channel) to 0.06 (a channel with stones and debris, 1/3 of vegetation)


if __name__ == "__main__":
    siz = (10,4)
    #
    te = np.zeros(siz)
    en = np.zeros((*siz,4))
    # water:
    wa = np.zeros(siz)
    wa[1,1] = 10
    #
    for t in range(1000):
        print("--- t=%d:"%t)
        print("wa:")
        print(wa)
        print(np.sum(wa))
        print("en:")
        print(en[:,:,0])
        print(en[:,:,1])
        print(en[:,:,2])
        print(en[:,:,3])
        print(np.sum(en))
        # Update
        wa,en = energy_transfer_step(wa,en,te,delta_t=0.001)
        input()
