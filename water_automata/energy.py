# nn = Manning coefficient of roughness - ranging from 0.01 (a clean and smooth channel) to 0.06 (a channel with stones and debris, 1/3 of vegetation)

def warp(z):
    w = np.zeros([z.shape[0]+2,z.shape[1]+2]+z.shape[2:])
    w[1:-1,1:-1] = z.copy()
    w[0,:] = w[-2,:]
    w[-1,:] = w[1,:]
    w[:,0] = w[:,-2]
    w[:,-1] = w[:,1]
    return w

def warp_delta(z1,z2):
    assert(z1.shape==z2.shape)
    d = np.zeros((*z.shape,4))
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
    w[0,:] += z[-1,;,2]
    # W
    w[:-1,:] += z[1:,:,3]
    w[-1,:] += z[0,:,3]
    #
    return w

def warp_transf_on_dir(z):
    wz = warp(z)
    t = np.zeros(z.shape)
    t[:,:,0] = wz[1:-1,:-2]
    t[:,:,1] = wz[1:-1,2:]
    t[:,:,2] = wz[:-2,1:-1]
    t[:,:,3] = wz[2:,1:-1]
    return t

def bounce(z):
    assert(z.shape[2]==4)
    bz = np.zeros(z.shape)
    bz[:,:,0] = z[:,:,1]
    bz[:,:,1] = z[:,:,0]
    bz[:,:,2] = z[:,:,3]
    bz[:,:,3] = z[:,:,2]
    return bz


def energy_transfer(wa,kine,te,gg=9.8,delta_t=1.0,hh=1.0,nn=0.01,pp=997,watransf=0.2):
    assert(kine.shape[:-1]==wa.shape)
    assert(te.shape==wa.shape)
    assert(kine.shape[2]==4) # WASD
    # Level and auxs
    wa_l = warp(wa)
    te_l = warp(te)
    le_l = wa_l+te_l
    wa4 = np.stack([wa,wa,wa,wa],axis=2)
    wa4_div = wa.copy()
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
    cw = np.stack([ca_d,ca_a,ca_s,ca_w],axis=2) - ca
    # Contact with terrain
    ct = wa4 - cw - ca
    # Contact with terrain coefficient:
    ct_coef = np.minimum(1,ct/wa4_div)
    # Discharge velocity
    divel = (1.0/nn)*wa4**(2/3.0)*(ca_d/hh)**0.5
    # Height transfer due discharge
    disch_w_t = delta_t*ca*divel/hh
    # Movement velocity
    mass = hh*hh*pp*wa4
    kivel = (2*kine/mass)**0.5
    # Height transfer due kinetic
    kine_w_t = delta_t*(ca+watransf*cw)*kivel/hh
    # Energy bounce:
    kine_e_b = (ct_coef)*kine*kivel/hh
    # Energy transfer:
    kine_e_t = (1.0-ct_coef)*kine*kivel/hh
    # @@@ HEIGHT TRANSFER @@@
    wa_t = np.minimum(disch_w_t+kine_w_t,wa4/5)
    n_wa = wa.copy()
    n_wa += warp_focus_sum(w_t) - np.sum(wa_t,axis=2)
    # Kinetic energy change due change in potential energy:
    kine_new = gg*(hh**2)*pp*warp_delta(wa,n_wa)*wa_t
    # New kinetic energy:
    n_kine = kine + warp_transf_on_dir(kine_new+kine_e_t+bounce(kine_e_b))
    n_kine -= np.sum(kine_e_t+kine_e_b,axis=2)
    n_kine = np.maximum(0,n_kine)
    #
    return wa,kine
