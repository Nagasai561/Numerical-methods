def rectangular_rule(function, a, b):
    return (function(a)*(b-a))
    
def midpoint_rule(function,a, b):
    return (b-a)*(function((a+b)/2))

def trapezoidal_rule(function, a, b):
    return (function(a) + function(b))*((b-a)/2)

def simpson_rule(function, a, b):
    return (function(a) + 4*function((a+b)/2)+ function(b))*((b-a)/6)

def composite_rules(function_rule, function, a, b, num_intervals):
    step_size = (b-a)/num_intervals
    result = 0
    for i in range(num_intervals):
        result += function_rule(function, a+i*step_size, a+i*(step_size+1))
    return result

def gaussian_rule(function, n):
    roots = [[-0.5773502692, 0.5773502692], [-0.7745966692, 0, 0.7745966692], [-0.8611363116, -0.3399810436, 0, 0.3399810436, 0.8611363116]]
    coefficients = [[1, 1], [0.55555556,0.88888889, 0.55555556], [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451]]
    root = roots[n-2]
    coefficient = coefficients[n-2]
    ans = 0
    for i in range(n):
        ans += coefficient[i]*function(root[i])
    return ans



