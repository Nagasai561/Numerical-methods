import root_finding

def euler_method(function, initial_x, initial_y, h, final_x):
    x = initial_x
    y = initial_y
    while(x < final_x):
        y += function(x, y)*(h)
        x += h
    return y

def taylor_method_order2(function, function_deri, initial_x, initial_y, h, final_x):
    x = initial_x
    y = initial_y
    while(x < final_x):
        y += h*(function(x, y)) + ((h*h)/2)*function_deri(x, y)
        x += h
    return y

def runge_kutta_method_order2(function, initial_x, initial_y, h, final_x, alpha, beta, a, b): #alpha and beta denote the fractions for choosing internal points
    x = initial_x                                                                               #'a' and 'b' denote the weights for the slopes
    y = initial_y
    while(x < final_x):
        k1 = h*function(x, y)
        k2 = h*function(x + alpha*h, y+beta*k1)
        y += a*k1 + b*k2
        x += h
    return y


def trapezoidal_rule(function, initial_x, initial_y, h, final_x, tolerance, max_iterations):
    x = initial_x
    y = initial_y
    while x < final_x:
        def f(y1):
            return y + (h/2)*(function(x, y) + function(x+h, y1))
        y = root_finding.fixed_point_method(f, y, tolerance, max_iterations)
        x += h
    return y


