import math as m

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

def composite_rule(rule, function, left, right, num_intervals):
    result = 0
    step = (right-left)/num_intervals
    for i in range(1, num_intervals+1):
        result += rule(function, left+(step*(i-1)), left+(step*i))
    return result

def gauss_quadrature_onepoint(function, left, right):
    def modified_input(input):
        return (input*(right-left) + (right+left))/2
    return 2*function(modified_input(left))

def gauss_quadrature_twopoints(function, left, right):
    def modified_input(input):
        return (input*(right-left) + (right+left))/2
    first_term = function(modified_input(0.577350269))
    second_term = function(modified_input(-0.577350269))
    return first_term + second_term
