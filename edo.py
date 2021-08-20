import numpy as np 


#Se definen los coeficientes que se van a utilizar en la integración.

b1 = 35/384                     ;           v1 = 5179/57600
b2 = 0                          ;           v2 = 0
b3 = 500/1113                   ;           v3 = 7571/16695
b4 = 125/192                    ;           v4 = 393/640
b5 = - 2187/6784                ;           v5 = - 92097/339200
b6 =  11/84                     ;           v6 = 187/2100
b7 = 0                          ;           v7 = 1/40


c2 , c3 , c4 , c5 , c6 , c7 = 1/5 , 3/10 , 4/5 , 8/9 , 1 , 1

a21 = 1/5
a31 = 3/40 ; a32 = 9/40
a41 = 44/45 ; a42 = - 56/15 ; a43 = 32/9
a51 = 19372/6561 ; a52 = - 25360/2187 ; a53 = 64448/6561 ; a54 = -212/729
a61=9017/3168 ; a62=-355/33 ; a63=46732/5247 ; a64=49/176 ; a65=-5103/18656
a71 = 35/384;a72=0;a73=500/1113;a74=125/192;a75=-2187/6784;a76=11/84



def rkck(f,x,y,h,*args):
    """Method that performes Runge-Kutta step of fifth order.
    Arguments:
        h : number, stepsize.
        f : callable as f(x,y) that returns a numpy array with the derivatives.
        x : point from wich to advance.
        y : state verctor in point x.
    Returns:
        tuple with the updated y value and the error.
    """
    k1 = h*f(x,y,*args)
    k2 = h*f(x+c2*h,y+a21*k1,*args)
    k3 = h*f(x+c3*h,y+a31*k1+a32*k2,*args)
    k4 = h*f(x+c4*h,y+a41*k1+a42*k2+a43*k3,*args)
    k5 = h*f(x+c5*h,y+a51*k1+a52*k2+a53*k3+a54*k4,*args)
    k6 = h*f(x+c6*h,y+a61*k1+a62*k2+a63*k3+a64*k4+a65*k5,*args)

    y1 = y + b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5 + b6*k6
    y2 = y + v1*k1 + v2*k2 + v3*k3 + v4*k4 + v5*k5 + v6*k6

    return y1 , abs(y2 - y1)

def rk(f,x,y,h,atol,rtol,hmax,hmin,*args):
    """Method that performes Runge-Kutta step of fifth order with error check.
    Arguments:
        f : callable as f(x,y) that returns a numpy array with the derivatives.
        x : point from wich to advance.
        y : state verctor in point x.
        h : number, attempted stepsize.
        atol : numpy array with desired absolute tolerance.
        rtol : numpy array with desired relative tolerance.
        hmax : largest allowed stepsize.
        hmin : shortest allowed stepsize.
    Returns:
        tuple containing the updated y , updated x , updated h in that order.
    """
    ynew , error = rkck(f,x,y,h,*args)
    error0 =  atol + rtol*abs(y)
    while any( error > error0 ):
        h = 0.8*h* np.sqrt(sum( (error0/error)**2 ) )**0.2
        if h >= hmax : 
            h = hmax 
            ynew , error = rkck(f,x,y,h,*args)  
            break
        if h <= hmin :
            h = hmin 
            ynew , error = rkck(f,x,y,h,*args)
            break
    return ynew , x + h , h















if __name__ == "__main__":
    from matplotlib.pylab import *

    x = array([1,0])

    def dxdt(t,x,w):
        return array([x[1]+x[0],-w**2* x[0]])

    v_val = [] ; x_val = [] ; t_val = []
    t = 0; h = 0.01
    while t <= 6*pi:
        x = rk(dxdt,t,x,h,array([0.001,0.001]),1e-3,0.1,0.01,2)[0]
        t_val.append(t)
        x_val.append(x[0])
        v_val.append(x[1])
        t += h
        print(t)

    figure()
    plot(t_val,x_val,label="x(t) [m]")
    plot(t_val,v_val,label="v(t) [m/s]")
    legend()
    xlabel("t[s]")
    show()

    figure()
    plot(x_val,v_val)
    xlabel("x(t) [m]") ; ylabel("v(t) [m/s]")
    show()
