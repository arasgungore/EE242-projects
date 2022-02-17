# Project 3: Secant and Bisection Methods



## Project Goals

In this project you will be implementing secant and bisection algorithms in order to solve 𝑓(𝑥) = 0 for
any given polynomial f.

Your program should take the coefficients of the function, initial guesses and the tolerance value as
command line arguments and return the resulting values of x as well as the numbers of iterations for
each method. For the bisection method, your program should give the warning that your algorithm
doesn’t work if the signs of initial guesses are same. You should implement both methods separately
first. Then you should use a hybrid method where you start with bisection methods for the first two
iterations and then continue with secant method for the rest of the iterations. Your program should print
out the number of iterations required for each of the 3 methods (i.e., bisection, secant, and hybrid). If
the number of iterations for each method exceeds 15, your program should be terminated and it should
give notes as “the number of iterations exceeded the threshold.”



## Project Details

- You have n+1 command line inputs for the coefficients of 𝑓(𝑥) = 𝑎<sub>𝑛</sub>𝑥<sup>𝑛</sup> + 𝑎<sub>𝑛-1</sub>𝑥<sup>𝑛-1</sup>
  + ⋯ + 𝑎 1 𝑥+𝑎 0 in the order from 𝑎𝑛 to 𝑎 0. Use dynamically allocated memory to store these.
- You have 3 more command line arguments for the initial guesses 𝑥<sub>0</sub>, 𝑥<sub>1</sub> (𝑥<sub>1</sub> > 𝑥<sub>0</sub>),
  and the tolerance value 𝑡𝑜𝑙.

**_Note:_** Please do not use functions or expressions for setting the output decimal precision. They truncate
or round the result. This may cause false positive or false negatives (i.e. Do not use setprecision
function). Show the result as it is.
