def line_fit(xs, ys): #takes in list of x and y coordinates and return a tuples (slope, y-intercept)
    xs_sum = 0
    ys_sum = 0
    xs_square_sum = 0
    xy_sum = 0
    n = len(xs)
    
    for x in xs:
        xs_sum += x
        xs_square_sum += x*x
    
    for y in ys:
        ys_sum += y

    for i in range(n):
        xy_sum += xs[i]*ys[i]

    m = (n*xy_sum-xs_sum*ys_sum)/(n*xs_square_sum-xs_sum*xs_sum)
    b = (ys_sum-m*xs_sum)/n

    return (m,b)

def binomial_coefficient(n, k): #calculates 'n choose k'.
    if(k < 0):
        return 1
    result = 1
    for i in range(0,k):
        result *= (n-i)
        result = result/(i+1)
    return result

