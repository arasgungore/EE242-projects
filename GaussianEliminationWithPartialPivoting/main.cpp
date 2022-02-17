// *** Header Imports ***
#include <iostream>			// Import cout.
//#include <iomanip>		// Import setf() and setprecision() functions.
#include <fstream>			// Import ifstream and ofstream.
#include <cmath>			// Import fabs() function.

#define E_MACH 1e-10		// Let's assume E_mach = 10^(-10).

// *** Function Prototypes ***
double* solveByGEPP(double **A, double *b, const short N);
void printCondNo(double **mat);
void deleteMatrix(double **mat, const short N);

// *** Main function ***
int main(int argc, char *argv[]) {
	try {
		short N = 0;							// N : Sizes of the A and b matrices, as A is a NxN square matrix and b is a Nx1 column matrix/vector.
		std::ifstream input_A(argv[1]);			// Input file stream declarations. (e.g. argv[1] = "A.txt", argv[2] = "b.txt")
		if(!input_A.is_open())					// Throw an exception if the input files with given filenames does not exist.
			throw argv[1];
		std::ifstream input_b(argv[2]);
		if(!input_b.is_open())
			throw argv[2];
		
		// Count the number of rows in the A matrix with a temporary local string variable "line".
		for(std::string line; std::getline(input_A, line); N++);
		
		double **A = new double*[N];			// Since we know N, now we can add the matrix declarations.
		double *b = new double[N];
		
		input_A.clear();						// Clears any fail bits from the input file stream such as the end-of-file bit.
		input_A.seekg(0, input_A.beg);			// Seeks to the beginning of the input file.
		
		for(short i = 0; i < N; i++) {
			A[i] = new double[N];				// Allocate memory for N integers to the i'th row of A matrix.
			for(short j = 0; j < N; j++)
				input_A >> A[i][j];				// Get each entry of the A matrix from the input file stream.
			input_b >> b[i];					// Get each entry of the b vector from the input file stream.
		}
		
		input_A.close();						// Close the input files.
		input_b.close();
		
		if(N == 2)								// If the size of A matrix is 2x2, then print out its condition numbers.
			printCondNo(A);
		double *x = solveByGEPP(A, b, N);		// Get the solution matrix.
		
		std::ofstream output_x("x.txt");		// Output file stream declaration.
		if(!output_x.is_open())					// Throw an exception if the output file with given filename cannot be generated.
			throw "x.txt";
		
	//	output_x.setf(std::ios::fixed, std::ios::floatfield);	// Fix floating point precision to 10 in the output file stream.
	//	output_x.setf(std::ios::showpoint);
	//	output_x.precision(10);
	//	std::cout.setf(std::ios::fixed, std::ios::floatfield);	// Fix floating point precision to 10 in the console output stream.
	//	std::cout.setf(std::ios::showpoint);
	//	std::cout.precision(10);
		
		std::cout << "x:" << std::endl;
		if(x == NULL) {							// If the x pointer is NULL, then the GEDD algorithm has failed to find a particular solution due to singularity of A matrix.
			std::cout << "\tNo particular solution.";
			output_x << "No particular solution.";
		}
		else {									// If the x pointer is not NULL, then the GEDD algorithm has successfully found a particular solution.
			for(short i = 0; i < N - 1; i++) {
				std::cout << "\t" << x[i] << std::endl;			// Print the solution both on the terminal and the output file.
				output_x << x[i] << std::endl;
			}
			std::cout << "\t" << x[N - 1];		// No newline at the end of the output.
			output_x << x[N - 1];
		}
		output_x.close();		// Close the output file.
		
		deleteMatrix(A, N);		// Free the memory of A matrix.
		delete [] b;			// Free the memory of b vector.
		return 0;				// Successful termination.
	}
	catch(const char *filename) {
		std::cerr << "[!ERROR!]: Could not open file \"" << filename << "\"." << std::endl;
		return 1;				// Unsuccessful termination.
	}
}


// *** Function Implementations ***

/**
 * Finds the solution to the linear system "A * x = b" using Gaussian Elimination with Partial Pivoting (GEPP) algorithm.
 * For more information about GEPP : http://www.math.iitb.ac.in/~neela/partialpivot.pdf
 * Time complexity: O(n^2)
 * @param A  Coefficient matrix A.
 * @param b  Column vector b.
 * @param N  Size of A and b matrices.
 * @return x Solution matrix x.
 *
 * procedure solveByGEPP(A, b, N)
 *     aug <- [A | b]
 *     for k = 0 to N - 1 step 1 do
 *         PivotSwap(aug, N, k)
 *         if max_pivot < E_mach then
 *             STOP
 *         end if
 *         for i = k + 1 to N - 1 step 1 do
 *             for j = k + 1 to N step 1 do
 *                 aug[i][j] <- aug[i][j] - aug[i][k] * aug[k][j] / aug[k][k]
 *             end for
 *             aug[i][k] <- 0
 *         end for
 *     end for
 *     for i = N - 1 to 0 step -1 do
 *         x[i] <- aug[i][N]
 *         for j = i + 1 to N - 1 step 1 do
 *             x[i] <- x[i] - aug[i][j] * x[j]
 *         end for
 *         x[i] <- x[i] / aug[i][i]
 *     end for
 *     return x
 * end procedure
 */
double* solveByGEPP(double **A, double *b, const short N) {
	double **aug = new double*[N];		// aug : Augmented matrix is a Nx(N+1) matrix obtained by appending the column vector b to A.
	double *x = NULL;					// x   : Nx1 solution matrix, which is equal to inv(A) * b.
	// To avoid compiler errors, I am going to use NULL instead of nullptr.
	
	for(short i = 0; i < N; i++) {		// Copy the entries of A to aug, then append vector b.
		aug[i] = new double[N + 1];
		for(short j = 0; j < N; j++)
			aug[i][j] = A[i][j];
		aug[i][N] = b[i];				// Aug[i] = [A_i1 A_i2 ... A_iN | b_i]
	}
	
	try {
		// Factorization.
		for(short k = 0; k < N; k++) {
			// In order to avoid general floating point arithmetic problems like underflow, overflow and roundoff errors...
			// ...we should choose the pivot as the entry that has the largest absolute value.
			short pivot_no = k;									// pivot_no  : Index of the row that contains max_pivot.
			double max_pivot = std::fabs(aug[k][k]);			// max_pivot : Absolute value of the pivot that has the largest absolute value.
			
			// Find the ideal pivot and its row number.
			for(short i = k + 1; i < N; i++) {		// Traverse column k to find the max_pivot.
				const double fabs_pivot_i = std::fabs(aug[i][k]);
				if(fabs_pivot_i > max_pivot) {
					max_pivot = fabs_pivot_i;
					pivot_no = i;
				}
			}
			
			// Throw an exception if the max_pivot is smaller than E_mach.
			// If the absolute value of the pivot is smaller than E_mach, the pivot is considered zero. Hence, the absence of a pivot implies singularity.
			if(max_pivot < E_MACH)
				throw max_pivot;
			
			// If the pivot isn't on the k'th row, swap rows pivot_no and k of the augmented matrix so that pivot is on the k'th row.
			if(pivot_no != k)
				for(short j = k; j < N + 1; j++)
					std::swap(aug[pivot_no][j], aug[k][j]);		// Swapping entries of these rows column by column.
			
			// After swapping rows we can move onto the next step, which is elimination of variables.
			for(short i = k + 1; i < N; i++) {
				const double coeff = aug[i][k] / aug[k][k];
				for(short j = k + 1; j < N + 1; j++)
					aug[i][j] -= coeff * aug[k][j];
				aug[i][k] = 0.0;
			}
		}
		
		x = new double[N];		// Allocate memory for N integers to the x vector.
		
		// Now that we have obtained upper triangular matrix U inside aug, we can solve the equation "U * x = b'" with back substitution. (aug = [U | b'] here)
		for(short i = N - 1; i >= 0; i--) {
			x[i] = aug[i][N];
			for(short j = i + 1; j < N; j++)
				x[i] -= aug[i][j] * x[j];
			x[i] /= aug[i][i];
		}
	}
	catch(const double &max_pivot) {
		// Print error messsages on the terminal.
		std::cerr << "[!ERROR!]: Absolute value of pivot is smaller than E_mach. (" << max_pivot << " < " << E_MACH << ")" << std::endl;
		std::cerr << "[!ERROR!]: Coefficient matrix A is not invertible. Solution x does not exist." << std::endl;
	}
	deleteMatrix(aug, N);		// Deallocate the augmented matrix.
	return x;
}

/**
 * Prints the condition numbers of 1 and infinity of the given 2x2 matrix on the terminal.
 * Time complexity: O(n^2)
 * @param mat Given matrix.
 *
 * pseudocode printCondNo(mat)
 *     N <- 2
 *     inv <- getInverse(mat, N)
 *     norm_1_mat <- 0, norm_1_inv <- 0, norm_inf_mat <- 0, norm_inf_inv <- 0
 *     for i = 0 to N step 1 do
 *         mat_col_sum <- 0, inv_col_sum <- 0
 *         for j = 0 to N step 1 do
 *             mat_col_sum <- abs(mat[j][i])
 *             inv_col_sum <- abs(inv[j][i])
 *         end for
 *         norm_1_mat <- max(norm_1_mat, mat_col_sum)
 *         norm_1_inv <- max(norm_1_inv, inv_col_sum)
 *     end for
 *     for i = 0 to N step 1 do
 *         mat_row_sum <- 0, inv_row_sum <- 0
 *         for j = 0 to N step 1 do
 *             mat_row_sum <- abs(mat[i][j])
 *             inv_row_sum <- abs(inv[i][j])
 *         end for
 *         norm_inf_mat <- max(norm_inf_mat, mat_row_sum)
 *         norm_inf_inv <- max(norm_inf_inv, inv_row_sum)
 *     end for
 *     cond_1 <- norm_1_mat * norm_1_inv
 *     cond_inf <- norm_inf_mat * norm_inf_inv
 *     print out cond_1, cond_inf
 * end procedure
 */
void printCondNo(double **mat) {
	const short N = 2;						// N   : Size of the given matrix, which is equal to 2.
	double **inv = new double*[N];			// inv : Inverse of the given matrix.
	for(short i = 0; i < N; i++)
		inv[i] = new double[N];
	
	const double det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];	// det : Determinant of the given matrix.
	inv[0][0] = mat[1][1];
	inv[0][1] = -mat[0][1];
	inv[1][0] = -mat[1][0];
	inv[1][1] = mat[0][0];
	for(short i = 0; i < N; i++)
		for(short j = 0; j < N; j++)
			inv[i][j] /= det;
	// For a matrix where N > 2, statements up to this point could be simplified to "double **inv = getInverse(mat, N);".
	// However, implementation of the getInverse function is out of the scope of the project, so I am going to leave a link to it here : https://www.geeksforgeeks.org/adjoint-inverse-matrix/
	
	double norm_1_mat = 0.0;		// norm_1_mat   : Vector 1-norm of the given matrix.
	double norm_1_inv = 0.0;		// norm_1_inv   : Vector 1-norm of the inverse of the given matrix.
	double norm_inf_mat = 0.0;		// norm_inf_mat : Vector infinity-norm of the given matrix.
	double norm_inf_inv = 0.0;		// norm_inf_inv : Vector infinity-norm of the inverse of the given matrix.
	for(short i = 0; i < N; i++) {
		double mat_col_sum = 0.0, inv_col_sum = 0.0;			// Vector 1-norm is equal to maximum absolute column sum.
		for(short j = 0; j < N; j++) {
			mat_col_sum += std::fabs(mat[j][i]);				// Calculate the absolute column sum of the given matrix for i'th column.
			inv_col_sum += std::fabs(inv[j][i]);				// Calculate the absolute column sum of the inverse of the given matrix for i'th column.
		}
		norm_1_mat = std::max(norm_1_mat, mat_col_sum);
		norm_1_inv = std::max(norm_1_inv, inv_col_sum);
	}
	for(short i = 0; i < N; i++) {
		double mat_row_sum = 0.0, inv_row_sum = 0.0;			// Vector infinity-norm is equal to maximum absolute row sum.
		for(short j = 0; j < N; j++) {
			mat_row_sum += std::fabs(mat[i][j]);				// Calculate the absolute row sum of the given matrix for i'th row.
			inv_row_sum += std::fabs(inv[i][j]);				// Calculate the absolute row sum of the inverse of the given matrix for i'th row.
		}
		norm_inf_mat = std::max(norm_inf_mat, mat_row_sum);
		norm_inf_inv = std::max(norm_inf_inv, inv_row_sum);
	}
	
	const double cond_1 = norm_1_mat * norm_1_inv;				// cond_1   : Condition number of 1 of the given matrix.
	const double cond_inf = norm_inf_mat * norm_inf_inv;		// cond_inf : Condition number of infinity of the given matrix.
	std::cout << "Condition number of 1: " << cond_1 << std::endl;
	std::cout << "Condition number of infinity: " << cond_inf << std::endl;
	deleteMatrix(inv, 2);			// Deallocate the inverse matrix.
}

/**
 * Frees/deallocates the memory of the given NxN matrix.
 * @param mat Given matrix.
 * @param N   Size of the given matrix.
 */
void deleteMatrix(double **mat, const short N) {
	for(short i = 0; i < N; i++)
		delete [] mat[i];			// Deallocate each row/pointer.
	delete [] mat;					// Deallocate the pointer to pointer.
}
