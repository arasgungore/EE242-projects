# Project 3: Secant and Bisection Methods

This project was assigned for the Numerical Methods for Electrical Engineering (EE 242) course in the Spring 2021 semester.



## Run on Terminal

<pre><code>g++ main.cpp -o test
test π<sub>n</sub> π<sub>π-1</sub> ... π<sub>0</sub> π₯<sub>0</sub> π₯<sub>1</sub> π‘ππ
(e.g. test 2 2 -7 1 -7 1.5 1.8 0.001)</code></pre>



## Project Goals

In this project, you will be implementing secant and bisection algorithms to solve π(π₯) = 0 for
any given polynomial f.

Your program should take the coefficients of the function, initial guesses, and the tolerance value as
command-line arguments and return the resulting values of x as well as the numbers of iterations for
each method. For the bisection method, your program should give the warning that your algorithm
doesnβt work if the signs of initial guesses are the same. You should implement both methods separately
first. Then you should use a hybrid method where you start with bisection methods for the first two
iterations and then continue with the secant method for the rest of the iterations. Your program should print
out the number of iterations required for each of the 3 methods (i.e., bisection, secant, and hybrid). If
the number of iterations for each method exceeds 15, your program should be terminated and it should
give notes as βthe number of iterations exceeded the threshold.β



## Project Details

- You have n+1 command-line inputs for the coefficients of π(π₯) = π<sub>π</sub>π₯<sup>π</sup> + π<sub>π-1</sub>π₯<sup>π-1</sup> + β― + π<sub>1</sub>π₯ + π<sub>0</sub>
in the order from π<sub>n</sub> to π<sub>0</sub>. Use dynamically allocated memory to store these.
- You have 3 more command line arguments for the initial guesses π₯<sub>0</sub>, π₯<sub>1</sub> (π₯<sub>1</sub> > π₯<sub>0</sub>),
  and the tolerance value π‘ππ.
