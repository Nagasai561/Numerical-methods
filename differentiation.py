def forward_difference(function, a, h):
    return (function(a+h)-function(a))/h

def backward_difference(function, a, h):
    return (function(a)-function(a-h))/h

def centered_difference(function, a, h):
    return (function(a+h)-function(a-h))/(2*h)