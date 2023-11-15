import root_finding
import linear_algebra
import math

REAL_COMPARISON_ERROR_TOLERANCE = 0.000001
FIXED_POINT_ITERATION_TOLERANCE = 0.0000001
MAX_FIXED_POINT_ITERATIONS = 100


#---------------------------------------Initial value problems---------------------------------



def euler_method(f, xy_o, h, final_x): #xy is [t, y] at initial point
    xy = xy_o.copy()
    while abs(xy[0]-final_x) > REAL_COMPARISON_ERROR_TOLERANCE:
        xy[1] += h*f(xy)
        xy[0] += h
    return xy[1]


def taylor_method_order2(f, ft_ffy, xy_o, h, final_x):
    xy = xy_o.copy()
    while abs(xy[0]-final_x) > REAL_COMPARISON_ERROR_TOLERANCE:
        xy[1] += h*f(xy)+h*(h/2)*ft_ffy(xy)
        xy[0] += h
    return xy[1]


def midpoint_method(f, xy_o, h, final_x): #runge_kutta_order2
    xy = xy_o.copy()
    while abs(xy[0]-final_x) > REAL_COMPARISON_ERROR_TOLERANCE:
        xy[1] += h*f([xy[0]+h/2, xy[1]+(h/2)*f(xy)])
        xy[0] += h
    return xy[1]


def modified_euler_method(f, xy_o, h, final_x):   #runge_kutta_order2
    xy = xy_o.copy()
    while abs(xy[0]-final_x) > REAL_COMPARISON_ERROR_TOLERANCE:
        xy[1] += (h/2)*f(xy)+(h/2)*f([xy[0]+h, xy[1]+h*f(xy)])
        xy[0] += h
    return xy[1]


def runge_kutta_order4(f, xy_o, h, final_x):
    xy = xy_o.copy()
    while abs(xy[0]-final_x) > REAL_COMPARISON_ERROR_TOLERANCE:
        k1 = h*f(xy)
        k2 = h*f([xy[0]+h/2, xy[1]+k1/2])
        k3 = h*f([xy[0]+h/2, xy[1]+k2/2])
        k4 = h*f([xy[0]+h, xy[1]+k3])
        xy[1] += (1/6)*(k1+2*k2+2*k3+k4)
        xy[0] += h
    return xy[1]

def adam_bashforth_order4(f, xy_o, h, final_x, *args): #if you are providing additional [ti, wi] 
    xy = xy_o.copy()                                    #you must pass it like [[t1, w1], [t2, w2], [t3, w3]]
    if(len(args) == 0):                                
        xy0 = xy
        xy1 = [xy[0]+h, runge_kutta_order4(f, xy, h, xy[0]+h)]
        xy2 = [xy[0]+2*h, runge_kutta_order4(f, xy, h, xy[0]+2*h)]
        xy3 = [xy[0]+3*h, runge_kutta_order4(f, xy, h, xy[0]+3*h)]
    else:
        xy0 = xy
        xy1 = args[0]
        xy2 = args[1]
        xy3 = args[2]

    while abs(xy3[0]-final_x) > REAL_COMPARISON_ERROR_TOLERANCE:
        temp = xy3.copy()
        xy3[1] += (h/24)*(55*f(xy3)-59*f(xy2)+37*f(xy1)-9*f(xy0))
        xy0 = xy1
        xy1 = xy2
        xy2 = temp
        xy3[0] += h
        

    return xy3[1]


def adam_moulton_order4(f, xy_o, h, final_x, *args):  #if you want to pass additional [ti, wi]
    xy = xy_o.copy()                                  #you must pass it like [[t1, w1], [t2, w2], [t3, w3]]
    if(len(args) == 0):                                
        xy0 = xy
        xy1 = [xy[0]+h, runge_kutta_order4(f, xy, h, xy[0]+h)]
        xy2 = [xy[0]+2*h, runge_kutta_order4(f, xy, h, xy[0]+2*h)]
    else:
        xy0 = xy
        xy1 = args[0]
        xy2 = args[1]

    while abs(xy2[0]-final_x) > REAL_COMPARISON_ERROR_TOLERANCE:
        temp = xy2.copy()
        def temp_fun(x):
            res = xy2[1]
            res += (h/24)*(9*f([xy2[0]+h, x])+19*f(xy2)-5*f(xy1)+f(xy0))
            return res
        xy2[1] = root_finding.fixed_point_method(temp_fun, xy2[1], FIXED_POINT_ITERATION_TOLERANCE, MAX_FIXED_POINT_ITERATIONS)
        xy2[0] += h
        xy0 = xy1
        xy1 = temp
    
    return xy2[1]


def predictor_corrector(f, xy_o, h, final_x, *args):     #if you want to pass additional [ti, wi]
    xy = xy_o.copy()                                    #you must pass it like [t1, w1], [t2, w2], [t3, w3]
    if(len(args) == 0):                                 
        xy0 = xy
        xy1 = [xy[0]+h, runge_kutta_order4(f, xy, h, xy[0]+h)]
        xy2 = [xy[0]+2*h, runge_kutta_order4(f, xy, h, xy[0]+2*h)]
        xy3 = [xy[0]+3*h, runge_kutta_order4(f, xy, h, xy[0]+3*h)]
    else:
        xy0 = xy
        xy1 = args[0]
        xy2 = args[1]
        xy3 = args[2]

    while abs(xy3[0]-final_x) > REAL_COMPARISON_ERROR_TOLERANCE:
        temp = xy3.copy()
        xy3[1] += (h/24)*(55*f(xy3)-59*f(xy2)+37*f(xy1)-9*f(xy0))
        xy0 = xy1
        xy1 = xy2
        xy2 = temp
        xy3[0] += h

    xy3[1] = xy2[1] + (h/24)*(9*f(xy3)+19*f(xy2)-5*f(xy1)+f(xy0))
    return xy3[1]

    
def runge_kutta_order4_for_multi_variable(f_list, xy_o, h, final_x): # xy_0 is [t0, x1, x2, ...] = u
    xy = xy_o.copy()                                    # fi is the function related to xi
    n = len(xy)                                          # each fi should take u as argument
    while abs(xy[0]-final_x) > REAL_COMPARISON_ERROR_TOLERANCE:
        temp = [0]*n

        k1 = [0]*n
        for i in range(1,n):
            k1[i] = h*(f_list[i-1](xy))
        
        k2 = [0]*n
        temp[0] = xy[0] + h/2
        for i in range(1, n):
            temp[i] = xy[i]+k1[i]/2
        for i in range(1, n):
            k2[i] = h*f_list[i-1](temp)
        
        k3 = [0]*n
        for i in range(1, n):
            k3[i] = xy[i] + k2[i]/2
        for i in range(1, n):
            k3[i] = h*f_list[i-1](temp)

        k4 = [0]*n
        temp[0] = xy[0] + h
        for i in range(1, n):
            temp[i] = xy[i]+k3[i]
        for i in range(1, n):
            k4[i] = h*f_list[i-1](temp)
        
        for i in range(1, n):
            xy[i] += (1/6)*(k1[i]+2*k2[i]+2*k3[i]+k4[i])

        xy[0] += h

    return xy


#y'' = p(x) y' + q(x) y + r(x)
# f should take [t, u1, u2] where u1 = y, u2 = y'
#xy_o = [t(0), y(0), y'(0)]
def runge_kutta_order4_for_second_order(f, xy_0, h, final_x):   
    y_dash = lambda u: u[2]
    f_list = [y_dash, f]
    return runge_kutta_order4_for_multi_variable(f_list, xy_0, h, final_x)


#--------------------------------------Boundary value problems-------------------------

#question to be solved y'' = p(x) y' + q(x) y + r(x) with y(a) = alpha and y(b) = beta
#ini_cond_1 = [a, alpha] and ini_cond_2 = [b, beta]
#func1 should be p(x) y' + q(x) y + r(x)
#func2 should be p(x) y' + q(x) y
#both should take [t, y, y']
def linear_shooting_method(func1, func2, ini_cond_1, ini_cond_2, h, final_x):
    sol1 = runge_kutta_order4_for_second_order(func1, [*ini_cond_1, 0], h, final_x)[1]
    b1 = runge_kutta_order4_for_second_order(func1,[*ini_cond_1, 0], h, ini_cond_2[0])[1]
    sol2 = runge_kutta_order4_for_second_order(func2, [ini_cond_1[0], 0, 1], h, final_x)[1]
    b2 = runge_kutta_order4_for_second_order(func2, [ini_cond_1[0], 0, 1], h, ini_cond_2[0])[1]
    res = sol1
    res += ((ini_cond_2[1]-b1)/b2)*sol2
    return res



#y'' = f(x, y, y')
#ini_cond_1 and ini_cond_2 are initial conditions respectively
#t_guess1 and t_guess2 are initial guesses
def non_linear_shooting_secant(f, ini_cond_1, ini_cond_2, h, final_x, t_guess1, t_guess2):
    def temp(t):
        return runge_kutta_order4_for_second_order(f, [*ini_cond_1, t], h, ini_cond_2[0])[1]-ini_cond_2[1]
    
    t = root_finding.secant_method(temp, t_guess1, t_guess2, FIXED_POINT_ITERATION_TOLERANCE, MAX_FIXED_POINT_ITERATIONS)
    return runge_kutta_order4_for_second_order(f, [*ini_cond_1, t], h, final_x)



#iterations is the number of times you want to improve your guess for t = y'(a)
#y'' = f(x, y, y')
#df_dy (a function) is partial differentiation of f wrt y: it takes in the parameters [x, y, y']
#df_dydash (a function) is partial differentiation of f wrt y': it takes in the paramenters [x, y, y']
def non_linear_shooting_newton(f, ini_cond_1, ini_cond_2, h, iterations, df_dy, df_dydash, t_guess, final_x):
    for iter in range(iterations):
        curr_sol = runge_kutta_order4_for_second_order(f, [*ini_cond_1, t_guess], h, ini_cond_2[0])
        def z(x):
            res = df_dy(curr_sol)*x[1]
            res += df_dydash(curr_sol)*x[2]
            return res
        z = runge_kutta_order4_for_second_order(z, [ini_cond_1[0], 0, 1], h, ini_cond_2[0])[1]
        t_guess = (t_guess)-(curr_sol[1]-ini_cond_2[1])/(z)

    return runge_kutta_order4_for_second_order(f, [*ini_cond_1, t_guess], h, final_x)


    


#y'' = p(x) y' + q(x) y + r(x)
#n is the size of matrix you want to create
#this divides the region [a, b] into n+1 equal intervals
def finite_difference_method(p, q, r, ini_cond_1, ini_cond_2, n):
    h = (ini_cond_2[0]-ini_cond_1[0])/(n+1)
    matrix = [[0]*n for i in range(n)]
    matrix[0][0] = 2+(h**2)*q(ini_cond_1[0]+h)
    matrix[0][1] = -1+(h/2)*(p(ini_cond_1[0]+h))
    matrix[n-1][n-1] = 2+(h**2)*q(ini_cond_2[0]-h)
    matrix[n-1][n-2] = -1-(h/2)*p(ini_cond_2[0]-h)
    for i in range(1, n-1):
        matrix[i][i-1] = -1-(h/2)*p(ini_cond_1[0]+(i+1)*h)
        matrix[i][i] = 2+(h**2)*(q(ini_cond_1[0]+(i+1)*h))
        matrix[i][i+1] = -1 + (h/2)*(p(ini_cond_1[0]+(i+1)*h))

    column_mat = [0]*n
    column_mat[0] = -(h**2)*(r(ini_cond_1[0]+h)) + (1 + (h/2)*(p(ini_cond_1[0]+h)))*ini_cond_1[1]
    column_mat[n-1] = -(h**2)*(r(ini_cond_2[0]-h)) + (1 - (h/2)*(p(ini_cond_2[0]-h)))*ini_cond_2[1]
    for i in range(1, n-1):
        column_mat[i] = -(h**2)*r(ini_cond_1[0]+(i+1)*h)

    return linear_algebra.gauss_elimination(matrix, column_mat)


    
    

