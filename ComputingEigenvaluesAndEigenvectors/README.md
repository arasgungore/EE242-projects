# Project 2: Computing Eigenvalues and Eigenvectors using Normalized
# Power Iteration Algorithm and Inverse Iteration Algorithm



## Project Goals

In this project, you will be implementing the normalized power iteration algorithm and
the inverse iteration algorithm to find the largest and the smallest eigenvalues and
corresponding eigenvectors of a given matrix A, where A is a real square matrix. Using
normalized power iteration, you can find the dominant eigenvalue of A with
the corresponding eigenvector. The inverse iteration algorithm is equivalent to the power
iteration method applied to 𝐴<sup>−1<sup>. Therefore, using inverse iteration algorithm the
smallest eigenvalue of A is the reciprocal of the dominant eigenvalue of 𝐴<sup>−1<sup> and the
corresponding eigenvector can be found. You should solve a system of linear equations
instead of finding 𝐴<sup>−1<sup> directly. (**Hint: You can use your first projects to solve the linear
system of equations, but this time you should write your program as an object-oriented one.**)

Your program should read A from an input file and output the dominant eigenvalue, its
corresponding eigenvector and the smallest eigenvalue with its corresponding
eigenvector as a text file.



## Programming Details

Your program should,

- have three command-line arguments for the parameters. (Command line arguments
can be thought of as the inputs of the main function.) The first argument is the name of
the file you read the matrix from, the second argument is the tolerance, which will be
used in the normalized power iteration algorithm and inverse iteration algorithm, and
the third argument is the name of your output file,
- use dynamically allocated memory to store the matrix,
- print out an error message and quit if it detects any problems,
- calculate the eigenvalues and corresponding eigenvector and write them in a text file.


## About Object-Oriented Programming

In this project, you will use matrices and matrix operations. In order to implement this
project as an object-oriented program; for example, you can declare a class (an object)
named Matrix, and implement all matrix operations such as multiplication, addition,
transpose, and solving a system of linear equations (LU factorization and backward
substitution) as methods of this class. This way, you don’t need to have a lot of complex
loops in your program.
