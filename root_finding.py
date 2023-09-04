def bisection_method(function, left_end, right_end, tolerance, max_iterations):
    iteration = 0
    while (iteration < max_iterations):
        iteration += 1
        midpoint = (left_end+right_end)//2
        value = function(midpoint)

        if(value == 0):
            return midpoint
        elif(value*function(left_end) > 0):
            left_end = midpoint
        else:
            right_end = midpoint

        if((right_end-left_end) < tolerance):
            return left_end

    return left_end

def fixed_point_method(function, initial_guess, tolerance, max_iterations):
    iteration = 0
    while(iteration < max_iterations):
        iteration += 1
        new_guess = function(initial_guess)

        if(abs(new_guess-initial_guess) < tolerance):
            return new_guess
        
        initial_guess = new_guess
    
    return initial_guess

def newton_method(function, function_derivative, initial_guess, tolerance, max_iterations):
    iteration = 0
    while(iteration < max_iterations):
        iteration += 1
        new_guess = initial_guess - (function(initial_guess)/function_derivative(initial_guess))
        
        if(abs(new_guess-initial_guess) < tolerance):
            return new_guess
        
        initial_guess = new_guess
    return initial_guess   
