from numpy import zeros, array , exp , pi , loadtxt
from lib.edo import rk
from time import process_time
from scipy.interpolate import InterpolatedUnivariateSpline as Spline

MS = 1.477

#Parametros para el Runge-Kutta de orden 5 con que se integrará.
rk_min_step = 0.01
rk_max_step = 0.01	 
rk_step = 0.002
g = None
pmin = None
atol = None 
rtol = None 
max_radius = 20



def derivatives(r,y,E,Pl,Rho,Nb):
    """Esta es la función que ofrecerá el miembro derecho en el sistema de ecuaciones diferenciales correspondientes a la métrica cilíndrica.
    Arguments:
        r : number, distance to the centre.
        y : numpy array containing the dependient variables.
        E : spline that interpolates energy density with the perpendicular pressure.
        Pl : spline that interpolates parallel pressure with the perpendicular pressure.
        Rho : spline that interpolates baryon density with the perpendicular pressure.
        Nb : spline that interpolates baryon concentration with the perpendicular pressure
    Returns:
        numpy array with the derivatives.
    """
    #y = (p,alpha,pfi,psi,pfi',psi',l')
    if r <= rk_min_step:
        return zeros(7,float)
    [p,alpha,pfi,psi,dpfi,dpsi,dl] = y
    pl = Pl(p) 
    e = E(p)
    exp2a = exp(2*alpha)
    dpdr = -p*(dpfi + dpsi) - dpfi*e + dpsi*pl   
    dalpha = dpsi + dpfi - 4*pi*r*exp2a*pl +  4*pi*r*exp2a*e
    ddpfi = 4*pi*exp2a*(e + pl + 2*p + r*dpfi*e) - dpfi*4*pi*r*pl*exp2a - dpfi/r
    ddpsi = -4*pi*exp2a*(e-r*dpsi*e+pl+r*dpsi*pl-2*p) - dpsi/r
    dl = 4*pi*r*exp2a*(e-2*p-pl)/MS   
    ret = array([dpdr,dalpha,dpfi,dpsi,ddpfi,ddpsi,dl])
    if any( abs(ret) > 1e+200) :
    	print("Overflow en derivadas en r = ",r," .")
    return ret


def observables(p0,pl0,e0,rho0,nb0,derivatives,E,Pl,Rho,Nb):
    """This functions computes the observables of a star according to its central parameters.
    Arguments:
        p0 : number, value of central perpendicular pressure.
        pl : number, value of central parallel pressure.
        e0 : number, value of central energy density.
        rho0 : number, value of central baryonic density.
        nb0 : number, value of central baryonic concentration.
        E : spline that interpolates energy density with the perpendicular pressure.
        Pl : spline that interpolates parallel pressure with the perpendicular pressure.
        Rho : spline that interpolates baryon density with the perpendicular pressure.
        Nb : spline that interpolates baryon concentration with the perpendicular pressure.
    Returns:
        numpy array containing central perpendicular pressure, parallel pressure, energy density, baryonic density, baryonic concentration, radius, mass/poar radius, central baryonic concentration.
    """
    global g , pmin
    y = array([p0,0,0,0,0,0,0])
    r = 0 ; dr = rk_step

    while not isinstance(y[0],complex):
        l = y[6]
        dy = derivatives(r,y,E,Pl,Rho,Nb)
        if y[0] + dy[0]*rk_max_step <= pmin or Pl(y[0]) <= 0 : break
        if r >= max_radius :
            print("Too many iterations")
            return None
        y , r , dr = rk(derivatives,r,y,dr,atol,rtol,rk_max_step,rk_min_step,E,Pl,Rho,Nb)
    
    return [p0/g , pl0/g,  e0 /g, rho0 , r , l, nb0 ]


def main(name_in,name_out,index_E,index_Pl,index_P,index_Rho,index_Nb,convertion_factor):
    global g , pmin , atol , rtol
    rtol = 1e-3
    if convertion_factor == 1:
        g  = 1.72e-13 #Este es el factor de conversión entre MeV⁴ de la unidades naturales
	      	#al km⁻² de las unidades geométricas.
    elif convertion_factor == 2:
        g = 1.32e-6# 1km^-2 = g * MeV/fm^3
    else:
        print("This program only works with these units.")

    start_time = process_time()
    dat = loadtxt(name_in)
    read_lines = len(dat)
    print("Read ",read_lines," lines.")

    E_dat = dat[:,index_E]*g
    Pl_dat = dat[:,index_Pl]*g
    P_dat = dat[:,index_P]*g
    Rho_dat = dat[:,index_Rho]
    Nb_dat = dat[:,index_Nb]

    #Se desechan aquellas lineas en que la presion perpendicular no es estrictamente creciente.
    init = 0
    while P_dat[init + 1] <= P_dat[init] and init + 1 < len(P_dat):
        init += 1
    if init != 0 : 
        print("Discarded first  ",init," lines.")
        init += 1

    E_dat = E_dat[init:]
    Pl_dat = Pl_dat[init:]
    P_dat = P_dat[init:]
    Rho_dat = Rho_dat[init:]
    Nb_dat = Nb_dat[init:]

    pmin = max(0,P_dat[0])
    atol = array([pmin,1e-3,1e-3,1e-3,1e-3,1e-6,1e-3])

    #Splines creation
    E = Spline(P_dat,E_dat)
    Pl = Spline(P_dat,Pl_dat)
    Rho = Spline(P_dat,Rho_dat)
    Nb = Spline(P_dat,Nb_dat)

    #Determination of first valid set of central parameters.
    init_val = 0
    for i in range(len(P_dat)):
        if P_dat[i] > 0 and Pl_dat[i] > 0:
            init_val = i
            break

    #File that will keep the observables.
    obs = open(name_out,"w")
    print("#p0  pl   e0  rho   r   m/z nb",file=obs)

    try:
        for i in range(init_val,len(P_dat)):
            if Rho_dat[i] < 1e+14 : continue
            print("Star # ",i)
            temp = observables(P_dat[i],Pl_dat[i],E_dat[i],Rho_dat[i],Nb_dat[i],derivatives,E,Pl,Rho,Nb)
            if temp :
                if temp[5] < 0 : break
                print(*temp,file=obs)

   
    
    finally :
        obs.close()

    time = process_time() - start_time
    


    return read_lines, init , time


if __name__ == "__main__":
#obtencion de datos y apertura de ficheros.
    name_in = input("Introducir nombre del fichero con datos a interpolar.\n")
    name_out = input("Introducir nombre del fichero de salida.\n")
    
    read_lines, descarded_lines , time = main(name_in,name_out,0,1,2,3,4,1)
    print("Read lines :",read_lines)
    print("Descarded lines :",descarded_lines)
    print("Took %d minuts and %.2lf seconds"%( time//60, time -60*time//60))
    input()


#E:\Facultad de Física\Segundo Año. Segundo Semestre\JCE\CalculoDeObservables\Python\lib\bosons_T0_B517_a1.dat
#E:\Facultad de Física\Segundo Año. Segundo Semestre\JCE\CalculoDeObservables\Python\lib\cbosons_T0_B517_a1.dat














