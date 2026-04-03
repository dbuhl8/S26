import numpy as np

# these are Finite difference schemes for fields taken from the simdat files
def FD4Z(field,nz,dz):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dz_field = np.zeros_like(field)
    dz_field[:,0,:,:] = (field[:,1,:,:]-field[:,0,:,:])/dz
    dz_field[:,1,:,:] = (field[:,2,:,:]-field[:,0,:,:])/(2*dz)
    dz_field[:,nz-2,:,:] = (field[:,nz-1,:,:]-field[:,nz-3,:,:])/(2*dz)
    dz_field[:,nz-1,:,:] = (field[:,nz-1,:,:]-field[:,nz-2,:,:])/dz
    # 4th order centered finite difference formula (inside)
    for i in range(nz-4):
        dz_field[:,i+2,:,:] = (1/(12*dz))*(-field[:,i+4,:,:]+8*field[:,i+3,:,:]\
                                           -8*field[:,i+1,:,:]+field[:,i,:,:])
    return dz_field

def FD6X(field,nx,dx):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dx_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    #dx_field[:,:,:,0] = (field[:,:,:,1]-field[:,:,:,0])/dx
    #dx_field[:,:,:,1] = (field[:,:,:,2]-field[:,:,:,0])/(2*dx)
    #dx_field[:,:,:,2] = (-field[:,:,:,4]+8*field[:,:,:,3] \
                        #-8*field[:,:,:,1]+field[:,:,:,0])/(12*dx)
    #dx_field[:,:,:,nx-3] = (-field[:,:,:,nx-1]+8*field[:,:,:,nx-2]\
                            #-8*field[:,:,:,nx-5]+field[:,:,:,nx-6])/(12*dx)
    #dx_field[:,:,:,nx-2] = (field[:,:,:,nx-1]-field[:,:,:,nx-3])/(2*dx)
    #dx_field[:,:,:,nx-1] = (field[:,:,:,nx-1]-field[:,:,:,nx-2])/dx
    dx_field[:,:,:,nx-1] = (1/(60*dx))*(-field[:,:,:,2]+9*field[:,:,:,1]\
                           -45*field[:,:,:,0]\
                           +45*field[:,:,:,nx-2]-9*field[:,:,:,nx-3]\
                           +field[:,:,:,nx-4])
    dx_field[:,:,:,nx-2] = (1/(60*dx))*(-field[:,:,:,1]+9*field[:,:,:,0]\
                           -45*field[:,:,:,nx-1]\
                           +45*field[:,:,:,nx-3]-9*field[:,:,:,nx-4]\
                           +field[:,:,:,nx-5])
    dx_field[:,:,:,nx-3] = (1/(60*dx))*(-field[:,:,:,0]+9*field[:,:,:,nx-1]\
                           -45*field[:,:,:,nx-2]\
                           +45*field[:,:,:,nx-4]-9*field[:,:,:,nx-5]\
                           +field[:,:,:,nx-6])

    dx_field[:,:,:,2] = (1/(60*dx))*(-field[:,:,:,5]+9*field[:,:,:,4]\
                           -45*field[:,:,:,3]\
                           +45*field[:,:,:,1]-9*field[:,:,:,0]\
                           +field[:,:,:,nx-1])
    dx_field[:,:,:,1] = (1/(60*dx))*(-field[:,:,:,4]+9*field[:,:,:,3]\
                           -45*field[:,:,:,2]\
                           +45*field[:,:,:,0]-9*field[:,:,:,nx-1]\
                           +field[:,:,:,nx-2])
    dx_field[:,:,:,0] = (1/(60*dx))*(-field[:,:,:,3]+9*field[:,:,:,2]\
                           -45*field[:,:,:,1]\
                           +45*field[:,:,:,nx-1]-9*field[:,:,:,nx-2]\
                           +field[:,:,:,nx-3])
    # 4th order centered finite difference formula (inside)
    for i in range(nx-6):
        dx_field[:,:,:,i+3] = (1/(60*dx))*(-field[:,:,:,i+6]+9*field[:,:,:,i+5]\
                                           -45*field[:,:,:,i+4]\
                                           +45*field[:,:,:,i+2]-9*field[:,:,:,i+1]\
                                           +field[:,:,:,i])
    return -dx_field

def FD6Y(field,ny,dy):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dy_field = np.zeros_like(field)
    #dy_field[:,:,0,:] = (field[:,:,1,:]-field[:,:,0,:])/dy
    #dy_field[:,:,1,:] = (field[:,:,2,:]-field[:,:,0,:])/(2*dy)
    #dy_field[:,:,2,:] = (-field[:,:,4,:]+8*field[:,:,3,:] \
                        #-8*field[:,:,1,:]+field[:,:,0,:])/(12*dy)
    #dy_field[:,:,ny-3,:] = (-field[:,:,ny-1,:]+8*field[:,:,ny-2,:]\
                            #-8*field[:,:,ny-5,:]+field[:,:,ny-6,:])/(12*dy)
    #dy_field[:,:,ny-2,:] = (field[:,:,ny-1,:]-field[:,:,ny-3,:])/(2*dy)
    #dy_field[:,:,ny-1,:] = (field[:,:,ny-1,:]-field[:,:,ny-2,:])/dy
    dy_field[:,:,ny-1,:] = (1/(60*dy))*(-field[:,:,2,:]+9*field[:,:,1,:]\
                           -45*field[:,:,0,:]\
                           +45*field[:,:,ny-2,:]-9*field[:,:,ny-3,:]\
                           +field[:,:,ny-4,:])
    dy_field[:,:,ny-2,:] = (1/(60*dy))*(-field[:,:,1,:]+9*field[:,:,0,:]\
                           -45*field[:,:,ny-1,:]\
                           +45*field[:,:,ny-3,:]-9*field[:,:,ny-4,:]\
                           +field[:,:,ny-5,:])
    dy_field[:,:,ny-2,:] = (1/(60*dy))*(-field[:,:,0,:]+9*field[:,:,ny-1,:]\
                           -45*field[:,:,ny-3,:]\
                           +45*field[:,:,ny-4,:]-9*field[:,:,ny-5,:]\
                           +field[:,:,ny-6,:])
    dy_field[:,:,2,:] = (1/(60*dy))*(-field[:,:,5,:]+9*field[:,:,4,:]\
                           -45*field[:,:,3,:]\
                           +45*field[:,:,1,:]-9*field[:,:,0,:]\
                           +field[:,:,ny-1,:])
    dy_field[:,:,1,:] = (1/(60*dy))*(-field[:,:,4,:]+9*field[:,:,3,:]\
                           -45*field[:,:,2,:]\
                           +45*field[:,:,0,:]-9*field[:,:,ny-1,:]\
                           +field[:,:,ny-2,:])
    dy_field[:,:,0,:] = (1/(60*dy))*(-field[:,:,3,:]+9*field[:,:,2,:]\
                           -45*field[:,:,1,:]\
                           +45*field[:,:,ny-1,:]-9*field[:,:,ny-2,:]\
                           +field[:,:,ny-3,:])
    # 4th order centered finite difference formula (inside)
    for i in range(ny-6):
        dy_field[:,:,i+3,:] = (1/(60*dy))*(-field[:,:,i+6,:]+9*field[:,:,i+5,:]\
                                           -45*field[:,:,i+4,:]\
                                           +45*field[:,:,i+2,:]-9*field[:,:,i+1,:]\
                                           +field[:,:,i,:])
    return -dy_field

def FD6Z(field,nz,dz):
    # takes in a field dataset with 4 dimensions (3 spatial)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dz_field = np.zeros_like(field)
    dz_field[:,nz-1,:,:] = (1/(60*dz))*(-field[:,2,:,:]+9*field[:,1,:,:]\
                           -45*field[:,0,:,:]\
                           +45*field[:,nz-2,:,:]-9*field[:,nz-3,:,:]\
                           +field[:,nz-4,:,:])
    dz_field[:,nz-2,:,:] = (1/(60*dz))*(-field[:,1,:,:]+9*field[:,0,:,:]\
                           -45*field[:,nz-1,:,:]\
                           +45*field[:,nz-3,:,:]-9*field[:,nz-4,:,:]\
                           +field[:,nz-5,:,:])
    dz_field[:,nz-3,:,:] = (1/(60*dz))*(-field[:,0,:,:]+9*field[:,nz-1,:,:]\
                           -45*field[:,nz-2,:,:]\
                           +45*field[:,nz-4,:,:]-9*field[:,nz-5,:,:]\
                           +field[:,nz-6,:,:])

    dz_field[:,2,:,:] = (1/(60*dz))*(-field[:,5,:,:]+9*field[:,4,:,:]\
                           -45*field[:,3,:,:]\
                           +45*field[:,1,:,:]-9*field[:,0,:,:]\
                           +field[:,nz-1,:,:])
    dz_field[:,1,:,:] = (1/(60*dz))*(-field[:,4,:,:]+9*field[:,3,:,:]\
                           -45*field[:,2,:,:]\
                           +45*field[:,0,:,:]-9*field[:,nz-1,:,:]\
                           +field[:,nz-2,:,:])
    dz_field[:,0,:,:] = (1/(60*dz))*(-field[:,3,:,:]+9*field[:,2,:,:]\
                           -45*field[:,1,:,:]\
                           +45*field[:,nz-1,:,:]-9*field[:,nz-2,:,:]\
                           +field[:,nz-3,:,:])
    # 4th order centered finite difference formula (inside)
    for i in range(nz-6):
        dz_field[:,i+3,:,:] = (1/(60*dz))*(-field[:,i+6,:,:]+9*field[:,i+5,:,:]\
                                           -45*field[:,i+4,:,:]\
                                           +45*field[:,i+2,:,:]-9*field[:,i+1,:,:]\
                                           +field[:,i,:,:])
    return dz_field

def iFD6Z(field,nz,dz):
    # takes in a field dataset with 4 dimensions (3 spatial)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dz_field = np.zeros_like(field)
    dz_field[nz-1,:,:] = (1/(60*dz))*(-field[2,:,:]+9*field[1,:,:]\
                           -45*field[0,:,:]\
                           +45*field[nz-2,:,:]-9*field[nz-3,:,:]\
                           +field[nz-4,:,:])
    dz_field[nz-2,:,:] = (1/(60*dz))*(-field[1,:,:]+9*field[0,:,:]\
                           -45*field[nz-1,:,:]\
                           +45*field[nz-3,:,:]-9*field[nz-4,:,:]\
                           +field[nz-5,:,:])
    dz_field[nz-3,:,:] = (1/(60*dz))*(-field[0,:,:]+9*field[nz-1,:,:]\
                           -45*field[nz-2,:,:]\
                           +45*field[nz-4,:,:]-9*field[nz-5,:,:]\
                           +field[nz-6,:,:])

    dz_field[2,:,:] = (1/(60*dz))*(-field[5,:,:]+9*field[4,:,:]\
                           -45*field[3,:,:]\
                           +45*field[1,:,:]-9*field[0,:,:]\
                           +field[nz-1,:,:])
    dz_field[1,:,:] = (1/(60*dz))*(-field[4,:,:]+9*field[3,:,:]\
                           -45*field[2,:,:]\
                           +45*field[0,:,:]-9*field[nz-1,:,:]\
                           +field[nz-2,:,:])
    dz_field[0,:,:] = (1/(60*dz))*(-field[3,:,:]+9*field[2,:,:]\
                           -45*field[1,:,:]\
                           +45*field[nz-1,:,:]-9*field[nz-2,:,:]\
                           +field[nz-3,:,:])
    # 4th order centered finite difference formula (inside)
    for i in range(nz-6):
        dz_field[i+3,:,:] = (1/(60*dz))*(-field[i+6,:,:]+9*field[i+5,:,:]\
                                           -45*field[i+4,:,:]\
                                           +45*field[i+2,:,:]-9*field[i+1,:,:]\
                                           +field[i,:,:])
    return dz_field

def iFD6X(field,nx,dx):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dx_field = np.zeros_like(field)
    dx_field[:,:,nx-1] = (1/(60*dx))*(-field[:,:,2]+9*field[:,:,1]\
                           -45*field[:,:,0]\
                           +45*field[:,:,nx-2]-9*field[:,:,nx-3]\
                           +field[:,:,nx-4])
    dx_field[:,:,nx-2] = (1/(60*dx))*(-field[:,:,1]+9*field[:,:,0]\
                           -45*field[:,:,nx-1]\
                           +45*field[:,:,nx-3]-9*field[:,:,nx-4]\
                           +field[:,:,nx-5])
    dx_field[:,:,nx-3] = (1/(60*dx))*(-field[:,:,0]+9*field[:,:,nx-1]\
                           -45*field[:,:,nx-2]\
                           +45*field[:,:,nx-4]-9*field[:,:,nx-5]\
                           +field[:,:,nx-6])

    dx_field[:,:,2] = (1/(60*dx))*(-field[:,:,5]+9*field[:,:,4]\
                           -45*field[:,:,3]\
                           +45*field[:,:,1]-9*field[:,:,0]\
                           +field[:,:,nx-1])
    dx_field[:,:,1] = (1/(60*dx))*(-field[:,:,4]+9*field[:,:,3]\
                           -45*field[:,:,2]\
                           +45*field[:,:,0]-9*field[:,:,nx-1]\
                           +field[:,:,nx-2])
    dx_field[:,:,0] = (1/(60*dx))*(-field[:,:,3]+9*field[:,:,2]\
                           -45*field[:,:,1]\
                           +45*field[:,:,nx-1]-9*field[:,:,nx-2]\
                           +field[:,:,nx-3])
    # 4th order centered finite difference formula (inside)
    for i in range(nx-6):
        dx_field[:,:,i+3] = (1/(60*dx))*(-field[:,:,i+6]+9*field[:,:,i+5]\
                                           -45*field[:,:,i+4]\
                                           +45*field[:,:,i+2]-9*field[:,:,i+1]\
                                           +field[:,:,i])
    return -dx_field

def iFD6Y(field,ny,dy):
    ## takes in a field dataset with 4 dimensions (3 spatial, 1 time)  # returns dataset (numpy) of same shape but computing the 1st derivative # with respect to x_com
    dy_field = np.zeros_like(field)
    #dy_field[:,:,0,:] = (field[:,:,1,:]-field[:,:,0,:])/dy
    #dy_field[:,:,1,:] = (field[:,:,2,:]-field[:,:,0,:])/(2*dy)
    #dy_field[:,:,2,:] = (-field[:,:,4,:]+8*field[:,:,3,:] \
                        #-8*field[:,:,1,:]+field[:,:,0,:])/(12*dy)
    #dy_field[:,:,ny-3,:] = (-field[:,:,ny-1,:]+8*field[:,:,ny-2,:]\
                            #-8*field[:,:,ny-5,:]+field[:,:,ny-6,:])/(12*dy)
    #dy_field[:,:,ny-2,:] = (field[:,:,ny-1,:]-field[:,:,ny-3,:])/(2*dy)
    #dy_field[:,:,ny-1,:] = (field[:,:,ny-1,:]-field[:,:,ny-2,:])/dy
    dy_field[:,ny-1,:] = (1/(60*dy))*(-field[:,2,:]+9*field[:,1,:]\
                           -45*field[:,0,:]\
                           +45*field[:,ny-2,:]-9*field[:,ny-3,:]\
                           +field[:,ny-4,:])
    dy_field[:,ny-2,:] = (1/(60*dy))*(-field[:,1,:]+9*field[:,0,:]\
                           -45*field[:,ny-1,:]\
                           +45*field[:,ny-3,:]-9*field[:,ny-4,:]\
                           +field[:,ny-5,:])
    dy_field[:,ny-2,:] = (1/(60*dy))*(-field[:,0,:]+9*field[:,ny-1,:]\
                           -45*field[:,ny-3,:]\
                           +45*field[:,ny-4,:]-9*field[:,ny-5,:]\
                           +field[:,ny-6,:])
    dy_field[:,2,:] = (1/(60*dy))*(-field[:,5,:]+9*field[:,4,:]\
                           -45*field[:,3,:]\
                           +45*field[:,1,:]-9*field[:,0,:]\
                           +field[:,ny-1,:])
    dy_field[:,1,:] = (1/(60*dy))*(-field[:,4,:]+9*field[:,3,:]\
                           -45*field[:,2,:]\
                           +45*field[:,0,:]-9*field[:,ny-1,:]\
                           +field[:,ny-2,:])
    dy_field[:,0,:] = (1/(60*dy))*(-field[:,3,:]+9*field[:,2,:]\
                           -45*field[:,1,:]\
                           +45*field[:,ny-1,:]-9*field[:,ny-2,:]\
                           +field[:,ny-3,:])
    # 4th order centered finite difference formula (inside)
    for i in range(ny-6):
        dy_field[:,i+3,:] = (1/(60*dy))*(-field[:,i+6,:]+9*field[:,i+5,:]\
                                           -45*field[:,i+4,:]\
                                           +45*field[:,i+2,:]-9*field[:,i+1,:]\
                                           +field[:,i,:])
    return -dy_field


def FD6X_xyslice(field,nx,dx):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dx_field = np.zeros_like(field)
    # foward / backwards finite difference formula (boundary)
    #dx_field[:,:,0] = (field[:,:,1]-field[:,:,0])/dx
    #dx_field[:,:,1] = (field[:,:,2]-field[:,:,0])/(2*dx)
    #dx_field[:,:,2] = (-field[:,:,4]+8*field[:,:,3] \
                      #-8*field[:,:,1]+field[:,:,0])/(12*dx)
    #dx_field[:,:,nx-3] = (-field[:,:,nx-1]+8*field[:,:,nx-2]\
                          #-8*field[:,:,nx-5]+field[:,:,nx-6])/(12*dx)
    #dx_field[:,:,nx-2] = (field[:,:,nx-1]-field[:,:,nx-3])/(2*dx)
    #dx_field[:,:,nx-1] = (field[:,:,nx-1]-field[:,:,nx-2])/dx


    dx_field[:,:,nx-1] = (1/(60*dx))*(-field[:,:,2]+9*field[:,:,1]\
                           -45*field[:,:,0]\
                           +45*field[:,:,nx-2]-9*field[:,:,nx-3]\
                           +field[:,:,nx-4])
    dx_field[:,:,nx-2] = (1/(60*dx))*(-field[:,:,1]+9*field[:,:,0]\
                           -45*field[:,:,nx-1]\
                           +45*field[:,:,nx-3]-9*field[:,:,nx-4]\
                           +field[:,:,nx-5])
    dx_field[:,:,nx-3] = (1/(60*dx))*(-field[:,:,0]+9*field[:,:,nx-1]\
                           -45*field[:,:,nx-2]\
                           +45*field[:,:,nx-4]-9*field[:,:,nx-5]\
                           +field[:,:,nx-6])

    dx_field[:,:,2] = (1/(60*dx))*(-field[:,:,5]+9*field[:,:,4]\
                           -45*field[:,:,3]\
                           +45*field[:,:,1]-9*field[:,:,0]\
                           +field[:,:,nx-1])
    dx_field[:,:,1] = (1/(60*dx))*(-field[:,:,4]+9*field[:,:,3]\
                           -45*field[:,:,2]\
                           +45*field[:,:,0]-9*field[:,:,nx-1]\
                           +field[:,:,nx-2])
    dx_field[:,:,0] = (1/(60*dx))*(-field[:,:,3]+9*field[:,:,2]\
                           -45*field[:,:,1]\
                           +45*field[:,:,nx-1]-9*field[:,:,nx-2]\
                           +field[:,:,nx-3])
    # 4th order centered finite difference formula (inside)
    for i in range(nx-6):
        dx_field[:,:,i+3] = (1/(60*dx))*(-field[:,:,i+6]+9*field[:,:,i+5]\
                                           -45*field[:,:,i+4]\
                                           +45*field[:,:,i+2]-9*field[:,:,i+1]\
                                           +field[:,:,i])
    return -dx_field

def FD6Y_xyslice(field,ny,dy):
    # takes in a field dataset with 4 dimensions (3 spatial, 1 time)
    # returns dataset (numpy) of same shape but computing the 1st derivative
    # with respect to x_comp
    dy_field = np.zeros_like(field)
    #dy_field[:,0,:] = (field[:,1,:]-field[:,0,:])/dy
    #dy_field[:,1,:] = (field[:,2,:]-field[:,0,:])/(2*dy)
    #dy_field[:,2,:] = (-field[:,4,:]+8*field[:,3,:] \
                        #-8*field[:,1,:]+field[:,0,:])/(12*dy)
    #dy_field[:,ny-3,:] = (-field[:,ny-1,:]+8*field[:,ny-2,:]\
                            #-8*field[:,ny-5,:]+field[:,ny-6,:])/(12*dy)
    #dy_field[:,ny-2,:] = (field[:,ny-1,:]-field[:,ny-3,:])/(2*dy)
    #dy_field[:,ny-1,:] = (field[:,ny-1,:]-field[:,ny-2,:])/dy
    dy_field[:,ny-1,:] = (1/(60*dy))*(-field[:,2,:]+9*field[:,1,:]\
                           -45*field[:,0,:]\
                           +45*field[:,ny-2,:]-9*field[:,ny-3,:]\
                           +field[:,ny-4,:])
    dy_field[:,ny-2,:] = (1/(60*dy))*(-field[:,1,:]+9*field[:,0,:]\
                           -45*field[:,ny-1,:]\
                           +45*field[:,ny-3,:]-9*field[:,ny-4,:]\
                           +field[:,ny-5,:])
    dy_field[:,ny-2,:] = (1/(60*dy))*(-field[:,0,:]+9*field[:,ny-1,:]\
                           -45*field[:,ny-3,:]\
                           +45*field[:,ny-4,:]-9*field[:,ny-5,:]\
                           +field[:,ny-6,:])
    dy_field[:,2,:] = (1/(60*dy))*(-field[:,5,:]+9*field[:,4,:]\
                           -45*field[:,3,:]\
                           +45*field[:,1,:]-9*field[:,0,:]\
                           +field[:,ny-1,:])
    dy_field[:,1,:] = (1/(60*dy))*(-field[:,4,:]+9*field[:,3,:]\
                           -45*field[:,2,:]\
                           +45*field[:,0,:]-9*field[:,ny-1,:]\
                           +field[:,ny-2,:])
    dy_field[:,0,:] = (1/(60*dy))*(-field[:,3,:]+9*field[:,2,:]\
                           -45*field[:,1,:]\
                           +45*field[:,ny-1,:]-9*field[:,ny-2,:]\
                           +field[:,ny-3,:])

    # 4th order centered finite difference formula (inside)
    for i in range(ny-6):
        dy_field[:,i+3,:] = (1/(60*dy))*(-field[:,i+6,:]+9*field[:,i+5,:]\
                                           -45*field[:,i+4,:]\
                                           +45*field[:,i+2,:]-9*field[:,i+1,:]\
                                           +field[:,i,:])
    return -dy_field

def rms(field):
    return np.sqrt((field**2).mean())
