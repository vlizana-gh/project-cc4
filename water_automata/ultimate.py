import numpy as np

def warp(z):
    assert(len(z.shape)==3)
    z_warp = np.zeros(z.shape)
    z_warp[0,:,1:]  = z[0,:,:-1]
    z_warp[0,:,0]   = z[0,:,-1]
    z_warp[1,:,:-1] = z[1,:,1:]
    z_warp[1,:,-1]  = z[1,:,0]
    z_warp[2,1:,:]  = z[2,:-1,:]
    z_warp[2,0,:]   = z[2,-1,:]
    z_warp[3,:-1,:] = z[3,1:,:]
    z_warp[3,-1,:]  = z[3,0,:]
    return z_warp

def inv(zw):
    assert(len(zw.shape)==3)
    return np.array([zw[1],zw[0],zw[3],zw[2]])

def ulti_step(w,u,v,d_t=0.01,d_x=1.0,gg=9.8,vv=1.0):
    # Discharges:
    w_4 = np.array([w,w,w,w])
    w_warp = warp(w_4)
    w_warp_inv = inv(w_warp)
    diff = np.maximum(0,w_4-w_warp_inv)
    mome = np.maximum(0,np.array([+u,-u,+v,-v]))
    disc = diff*vv*d_t/d_x**2 + mome*d_t/d_x
    # NOTE: discharge annulation! should prevent invalid roots.
    disc[:2] -= np.min(disc[:2],axis=0)
    disc[2:] -= np.min(disc[2:],axis=0)
    disc_warp = warp(disc)
    # Speeds:
    w_div = np.copy(w)
    w_div[w_div==0] = np.inf
    s_x = u/w_div
    s_y = v/w_div
    s_x_4 = np.array([s_x,s_x,s_x,s_x])
    s_y_4 = np.array([s_y,s_y,s_y,s_y])
    # Total incoming discharge
    tw_in = np.sum(disc_warp,axis=0)
    tw_in_4 = np.array([tw_in,tw_in,tw_in,tw_in])
    # Total outgoing discharge
    tw_out = np.sum(disc,axis=0)
    tw_out_4 = np.array([tw_out,tw_out,tw_out,tw_out])
    tw_out_warp = warp(tw_out_4)
    # Incoming mean height change
    k_inc = (w_warp-w_4)-0.5*tw_out_warp+tw_out_4-0.5*tw_in_4
    # Momentum transference
    u_tran = disc*s_x_4
    v_tran = disc*s_y_4
    # Momentum generation:
    u_gen = np.zeros(w.shape)
    u_gen += disc_warp[0]*((np.maximum(0,2*gg*k_inc[0]+s_x**2))**0.5-s_x)
    u_gen += disc_warp[1]*(-(np.maximum(0,2*gg*k_inc[1]+s_x**2))**0.5-s_x)
    #
    v_gen = np.zeros(w.shape)
    v_gen += disc_warp[2]*((np.maximum(0,2*gg*k_inc[2]+s_y**2))**0.5-s_y)
    v_gen += disc_warp[3]*(-(np.maximum(0,2*gg*k_inc[3]+s_y**2))**0.5-s_y)
    # Water update:
    neo_w = w+tw_in-tw_out
    # Momentum update:
    neo_u = u+u_gen+np.sum(warp(u_tran)-u_tran,axis=0)
    neo_v = v+v_gen+np.sum(warp(v_tran)-v_tran,axis=0)
    return neo_w,neo_u,neo_v


# TEST:
if __name__ == "__main__":
    siz = (5,3)
    u = np.zeros(siz)
    v = np.zeros(siz)
    w = np.zeros(siz)
    w[1:3,1:3] = 5
    #
    for t in range(1000):
        print("--- t=%d:"%t)
        print("w:")
        print(w)
        print(np.sum(w))
        print("u:")
        print(u)
        print("v:")
        print(v)
        # Update
        w,u,v = ulti_step(w,u,v)
        input()
