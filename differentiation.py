def forward_difference(function, a, h):
    return (function(a+h)-function(a))/h

def backward_difference(function, a, h):
    return (function(a)-function(a-h))/h

def centered_difference(function, a, h):
    return (function(a+h)-function(a-h))/(2*h)

#three point formulas
def three_point_endpoint(function, a, h):
    num = -3*function(a)+4*function(a+h)-function(a+2*h)
    return num/(2*h)

def three_point_midpoint(function, a, h):
    num = function(a+h)-function(a-h)
    return num/(2*h)


#five point formulas
def five_point_endpoint(function, a, h):
    num = -25*function(a)+48*function(a+h)-36*function(a+2*h)+16*function(a+3*h)-3*function(a+4*h)
    return num/(12*h)

def five_point_midpoint(function, a, h):
    num = function(a-2*h)-8*function(a-h)+8*function(a+h)-function(a+2*h)
    return num/(12*h)


#second derivative midpoint formula
def second_der_midpoint(function, a, h):
    num = function(a-h)-2*function(a)+function(a+h)
    return num/(h**2)
