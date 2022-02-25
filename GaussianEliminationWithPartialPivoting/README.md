# Project 1: Gaussian Elimination with Partial Pivoting


## Project Goals

In this project, you will be implementing the Gaussian elimination algorithm with partial
pivoting together with backward substitution to solve Ax = b, where A is an n by n square
matrix.

Your program should read A and b from two input text files and output the solution x as a text
file.



## Programming Details

Your program should,

- read A matrix from A.txt file and b vector from b.txt file. Each line in a file represents a row,
- use dynamically allocated memory to store the matrix and the vector,
- print out an error message and quit if it detects that A is singular, (Don’t forget to consider the
machine precision while detecting singularity.)
- print out the elements of the solution x (in the correct order) and write them in a text file.



## The Case of High Condition Numbers

- In order not to add a great amount of burden onto this project, you do not have to find the
condition numbers for every given matrix. However, in the case of 2x2 matrices, your program
is required to output the condition numbers using 1 - norm and infinity-norm. (You can print
these on the screen; they are not required in a text file.)



## Run on Terminal

```sh
g++ main.cpp -o test
test A.txt b.txt
```
