def rectangular_rule(function, left, right):
    value = function(left)
    h = (right-left)
    return value*h

def midpoint_rule(function, left, right):
    middle = (left+right)/2
    return (right-left)*function(middle)


def trapezoidal_rule(function, left, right):
    sum_heights = function(left)+function(right)
    h = (right-left)
    return (h/2)*sum_heights

def simpson_rule(function, left, right):
    middle = (left+right)/2
    h = (right-left)/2
    expr = function(left) + 4*function(middle) + function(right)
    return expr*(h/3)

def composite_trapezoidal_rule(function, left, right, num_intervals):
    result = 0
    step = (right-left)/num_intervals
    for i in range(1, num_intervals+1):
        result += trapezoidal_rule(function, left+(step*(i-1)), left+(step*(i)))
    return result

def composite_simpson_rule(function, left, right, num_intervals):
    result = 0
    step = (right-left)/num_intervals
    for i in range(1, num_intervals+1):
        result += simpson_rule(function, left+(step*(i-1)), left+(step*(i)))
    return result
