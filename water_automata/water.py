import numpy as np

def momentum_transfer_step(wa,mo_x,mo_y,te,alpha,gamma,delta_t=1.0,h=1.0,density=1.0):
    # Level
    le = wa+te
    # Copy of water used to avoid 0/0 errors.
    wa_div = wa.copy()
    wa_div[wa_div==0] = 1
    # Column difference
    diff_d = np.maximum(0,np.minimum(le[:,:-1]-le[:,1:],wa[:,:-1]))
    diff_a = np.maximum(0,np.minimum(le[:,1:]-le[:,:-1],wa[:,1:]))
    diff_s = np.maximum(0,np.minimum(le[:-1,:]-le[1:,:],wa[:-1,:]))
    diff_w = np.maximum(0,np.minimum(le[1:,:]-le[:-1,:],wa[1:,:]))
    # Column difference coefficent
    diff_coef_d = np.minimum(1.0,diff_d/wa_div[:,:-1])
    diff_coef_a = np.minimum(1.0,diff_a/wa_div[:,1:])
    diff_coef_s = np.minimum(1.0,diff_s/wa_div[:-1,:])
    diff_coef_w = np.minimum(1.0,diff_w/wa_div[1:,:])
    # Column contact
    cont_d = np.maximum(0,np.minimum(le[:,:-1]-te[:,1:],wa[:,:-1]))
    cont_a = np.maximum(0,np.minimum(le[:,1:]-te[:,:-1],wa[:,1:]))
    cont_s = np.maximum(0,np.minimum(le[:-1,:]-te[1:,:],wa[:-1,:]))
    cont_w = np.maximum(0,np.minimum(le[1:,:]-te[:-1,:],wa[1:,:]))
    # Discharge
    disch_d = alpha*delta_t*diff_d
    disch_a = alpha*delta_t*diff_a
    disch_s = alpha*delta_t*diff_s
    disch_w = alpha*delta_t*diff_w
    # Flow
    flow_d = diff_coef_d*(delta_t*np.maximum(0,mo_x[:,:-1])/h**3/density)
    flow_a = diff_coef_a*(delta_t*np.maximum(0,-mo_x[:,1:])/h**3/density)
    flow_s = diff_coef_s*(delta_t*np.maximum(0,mo_y[:-1,:])/h**3/density)
    flow_w = diff_coef_w*(delta_t*np.maximum(0,-mo_y[1:,:])/h**3/density)
    # Final transference
    wa_tran_d = np.minimum(disch_d+flow_d,wa[:,:-1]/5.0)
    wa_tran_a = np.minimum(disch_a+flow_a,wa[:,1:]/5.0)
    wa_tran_s = np.minimum(disch_s+flow_s,wa[:-1,:]/5.0)
    wa_tran_w = np.minimum(disch_w+flow_w,wa[1:,:]/5.0)
    # Transference coefficent
    wa_tran_coef_d = np.minimum(1.0,wa_tran_d/wa_div[:,:-1])
    wa_tran_coef_a = np.minimum(1.0,wa_tran_a/wa_div[:,1:])
    wa_tran_coef_s = np.minimum(1.0,wa_tran_s/wa_div[:-1,:])
    wa_tran_coef_w = np.minimum(1.0,wa_tran_w/wa_div[1:,:])
    # Momentum transference
    mo_x_tran_d = mo_x[:,:-1]*wa_tran_coef_d
    mo_x_tran_a = mo_x[:,1:]*wa_tran_coef_a
    mo_x_tran_s = mo_x[:-1,:]*wa_tran_coef_s
    mo_x_tran_w = mo_x[1:,:]*wa_tran_coef_w
    #
    mo_y_tran_d = mo_y[:,:-1]*wa_tran_coef_d
    mo_y_tran_a = mo_y[:,1:]*wa_tran_coef_a
    mo_y_tran_s = mo_y[:-1,:]*wa_tran_coef_s
    mo_y_tran_w = mo_y[1:,:]*wa_tran_coef_w
    # Update momentums due transference
    mo_x_f = mo_x.copy()
    mo_y_f = mo_y.copy()
    #
    mo_x_f[:,:-1] += mo_x_tran_a-mo_x_tran_d
    mo_x_f[:,1:] += mo_x_tran_d-mo_x_tran_a
    mo_x_f[:-1,:] += mo_x_tran_w-mo_x_tran_s
    mo_x_f[1:,:] += mo_x_tran_s-mo_x_tran_w
    #
    mo_y_f[:,:-1] += mo_y_tran_a-mo_y_tran_d
    mo_y_f[:,1:] += mo_y_tran_d-mo_y_tran_a
    mo_y_f[:-1,:] += mo_y_tran_w-mo_y_tran_s
    mo_y_f[1:,:] += mo_y_tran_s-mo_y_tran_w
    # Update water levels
    wa_f = wa.copy()
    wa_f[:,:-1] += wa_tran_a-wa_tran_d
    wa_f[:,1:] += wa_tran_d-wa_tran_a
    wa_f[:-1,:] += wa_tran_w-wa_tran_s
    wa_f[1:,:] += wa_tran_s-wa_tran_w
    # TODO: Generated momentum
    # Return
    return wa_f, mo_x_f, mo_y_f

if __name__ == "__main__":
    siz = (5,4)
    alpha = 0.1
    gamma = 0.3
    #
    te = np.zeros(siz)
    mo_x = np.zeros(siz)
    mo_y = np.zeros(siz)
    # water:
    wa = np.zeros(siz)
    wa[1,1] = 100
    mo_x[1,1] = 10
    for t in range(1000):
        print("t=%d:"%t)
        print(wa)
        print(mo_x)
        print(mo_y)
        # Update
        wa,mo_x,mo_y = momentum_transfer_step(wa,mo_x,mo_y,te,alpha,gamma)
        raw_input()
