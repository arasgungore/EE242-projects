// *** Header Imports ***
#include <iostream>			// Import cout.
//#include <iomanip>		// Import setf() and setprecision() functions.
#include <fstream>			// Import ifstream and ofstream.
#include <cmath>			// Import fabs() function.
#include <cstdlib>			// Import atof() function.
#include <random>			// Import uniform_real_distribution and default_random_engine.

#define ITERATION_MAX 10000			// Max number of iterations allowed.

// *** Class Declaration ***
class Matrix;

// Matrix class.
class Matrix {
private:
	int row, col;
	double **mat;
public:
	Matrix(const int row, const int col = 1);
	Matrix(const int row, const int col, double **mat);
	Matrix(const Matrix &M);
	~Matrix();
	bool hasSameDimensions(const Matrix &M) const;
	double getEntry(const int r, const int c = 0) const;
	Matrix operator + (const Matrix &M) const;
	Matrix operator - (const Matrix &M) const;
	Matrix operator * (const Matrix &M) const;
	Matrix operator * (const double scalar) const;
	Matrix operator / (const double scalar) const;
	Matrix& operator = (const Matrix &M);
	Matrix getTranspose() const;
	double getDotProduct(const Matrix &M) const;
	double getInfinityNorm() const;
	void printMatrix(const char* name) const;
	Matrix getDeflatedMatrix(const Matrix eigenvector, const double tol) const;
	std::pair<double, Matrix> normPowIteration(const double tol) const;
};


// *** Class Implementation ***

/**
 * Matrix constructor with 2 parameters.
 */
Matrix::Matrix(const int row, const int col) : row(row), col(col), mat(new double*[row]) {
	for(int i = 0; i < row; i++)
		mat[i] = new double[col]();
}

/**
 * Matrix constructor with 3 parameters.
 */
Matrix::Matrix(const int row, const int col, double **mat) : row(row), col(col), mat(mat) { }

/**
 * Matrix copy constructor.
 */
Matrix::Matrix(const Matrix &M) : row(M.row), col(M.col), mat(new double*[M.row]) {
	for(int i = 0; i < row; i++) {
		mat[i] = new double[col];
		for(int j = 0; j < col; j++)
			mat[i][j] = M.mat[i][j];
	}
}

/**
 * Matrix destructor.
 */
Matrix::~Matrix() {
	for(int i = 0; i < row; i++)
		delete [] mat[i];			// Deallocate each row/pointer.
	delete [] mat;					// Deallocate the pointer to pointer.
}

/**
 * Checks if the two matrices have the same dimensions or not.
 */
bool Matrix::hasSameDimensions(const Matrix &M) const {
	return row == M.row && col == M.col;
}

/**
 * Returns the entry at row "r" and column "c".
 */
double Matrix::getEntry(const int r, const int c) const {
	try {
		if(r < 0 || r >= row)
			throw r;
		if(c < 0 || c >= col)
			throw c;
	} catch(const int &dim) {
		std::cerr << "[!ERROR!]: Dimension " << dim << " is out of boundaries." << std::endl;
	}
	
	return mat[r][c];
}

/**
 * Matrix addition operator. Adds the given matrix to this matrix and returns the result.
 */
Matrix Matrix::operator + (const Matrix &M) const {
	try {
		if(!hasSameDimensions(M))
			throw;
	} catch(...) {
		std::cerr << "[!ERROR!]: Matrices have different dimensions. Cannot do matrix addition." << std::endl;
	}
	
	Matrix sum(*this);
	for(int i = 0; i < row; i++)
		for(int j = 0; j < col; j++)
			sum.mat[i][j] += M.mat[i][j];
	return sum;
}

/**
 * Matrix substraction operator. Substracts the given matrix from this matrix and returns the result.
 */
Matrix Matrix::operator - (const Matrix &M) const {
	try {
		if(!hasSameDimensions(M))
			throw;
	} catch(...) {
		std::cerr << "[!ERROR!]: Matrices have different dimensions. Cannot do matrix subtraction." << std::endl;
	}
	
	Matrix diff(*this);
	for(int i = 0; i < row; i++)
		for(int j = 0; j < col; j++)
			diff.mat[i][j] -= M.mat[i][j];
	return diff;
}

/**
 * Matrix multiplication operator. Multiplies this matrix with the given matrix and returns the result.
 */
Matrix Matrix::operator * (const Matrix &M) const {
	try {
		if(col != M.row)
			throw;
	} catch(...) {
		std::cerr << "[!ERROR!]: Matrix dimensions do not match. Cannot do matrix multiplication." << std::endl;
	}
	
	Matrix prod(row, M.col);
	for(int i = 0; i < row; i++)
		for(int j = 0; j < M.col; j++)
			for(int k = 0; k < col; k++)
				prod.mat[i][j] += mat[i][k] * M.mat[k][j];
	return prod;
}

/**
 * Scalar multiplication operator. Multiplies the matrix with the given scalar number and returns the result.
 */
Matrix Matrix::operator * (const double scalar) const {
	Matrix prod(*this);
	for(int i = 0; i < row; i++)
		for(int j = 0; j < col; j++)
			prod.mat[i][j] *= scalar;
	return prod;
}

/**
 * Scalar division operator. Divides the matrix by the given scalar number and returns the result.
 */
Matrix Matrix::operator / (const double scalar) const {
	Matrix quot(*this);
	for(int i = 0; i < row; i++)
		for(int j = 0; j < col; j++)
			quot.mat[i][j] /= scalar;
	return quot;
}

/**
 * Assignment operator. Assigns the given matrix to this matrix.
 */
Matrix& Matrix::operator = (const Matrix &M) {
	for(int i = 0; i < row; i++)
		delete [] mat[i];
	delete [] mat;
	
	row = M.row;
	col = M.col;
	mat = new double*[row];
	
	for(int i = 0; i < row; i++) {
		mat[i] = new double[col];
		for(int j = 0; j < col; j++)
			mat[i][j] = M.mat[i][j];
	}
	return *this;
}

/**
 * Returns the transpose of the matrix.
 */
Matrix Matrix::getTranspose() const {
	Matrix transpose(col, row);
	for(int i = 0; i < row; i++)
		for(int j = 0; j < col; j++)
			transpose.mat[j][i] = mat[i][j];
	return transpose;
}

/**
 * Returns the dot product of the two vectors (Nx1 matrices).
 */
double Matrix::getDotProduct(const Matrix &M) const {
	if(col != 1 || M.col != 1)
		exit(1);
	double sum = 0.0;
	for(int i = 0; i < row; i++)
		sum += mat[i][0] * M.mat[i][0];
	return sum;
}

/**
 * Returns the infinity norm of the matrix.
 */
double Matrix::getInfinityNorm() const {
	double max_absolute_row_sum = 0.0;
	for(int i = 0; i < row; i++) {
		double absolute_row_sum = 0.0;
		
		for(int j = 0; j < col; j++)
			absolute_row_sum += std::fabs(mat[i][j]);
		
		max_absolute_row_sum = std::max(max_absolute_row_sum, absolute_row_sum);
	}
	return max_absolute_row_sum;
}

/**
 * Prints the matrix on the terminal.
 */
void Matrix::printMatrix(const char* name) const {
	std::cout << name << ":" << std::endl;
	for(int i = 0; i < row; i++) {
		std::cout << "\t[";
		for(int j = 0; j < col; j++)
			std::cout << " " << mat[i][j];
		std::cout << " ]" << std::endl;
	}
}

/**
 * Returns the deflated matrix.
 *
 * Pseudocode:
 *
 * procedure getDeflatedMatrix(A, u)
 *     local p, A_p
 *     # compute p value
 *     p <- p_index(u)
 *     # pth row of A
 *     A_p <- A[p]
 *     # output
 *     return A - (u . ap) / u[p]
 * end procedure
 */
Matrix Matrix::getDeflatedMatrix(const Matrix eigenvector, const double tol) const {
	int p = 0;
	for(; std::fabs(eigenvector.mat[p][0]) <= tol; p++);		// p corresponds to the entry which has an above tolerance value.
	
	Matrix pth_row(1, col);
	for(int i = 0; i < col; i++)		// Get the p'th row of this matrix.
		pth_row.mat[0][i] = mat[p][i];
	
	Matrix deflated_matrix = *this - (eigenvector * pth_row) / eigenvector.mat[p][0];

	return deflated_matrix;
}

/**
 * Uses normalized power iteration algorithm to calculate the dominant eigenvalue and its corresponding eigenvector.
 * 
 * Pseudocode:
 *
 * procedure normPowIteration(A, tol)
 *     x0 <- arbitrary nonzero vector
 *     for k = 1, 2, ...
 *         # generate next vector
 *         y_k <- A * x_(k-1)
 *         # normalize
 *         x_k <- y_k / ||y_k||_inf
 *     end for
 * end procedure
 */
std::pair<double, Matrix> Matrix::normPowIteration(const double tol) const {
	Matrix v(row);
	
	// Generate a random vector to start with.
	std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
	std::default_random_engine rand_eng;
	for(int i = 0; i < row; i++)
		v.mat[i][0] = uniform_dist(rand_eng);

	Matrix w = *this * v;		// First iteration.
	double eigenvalue = w.getInfinityNorm();
	v = w / eigenvalue;
	Matrix prev_v = v;
	double prev_eigenvalue = eigenvalue;

	for(int i = 1; i < ITERATION_MAX; i++) {	// Other iterations are performed until tolerance is broken or the iteration limit is reached.
		w = *this * v;
		eigenvalue = w.getInfinityNorm();
		v = w / eigenvalue;
		
		// Break if the the new and old eigenvalues found in the last 2 iterations are close enough so that their absolute difference is below the tolerance level.
		if(std::fabs(eigenvalue - prev_eigenvalue) < tol)
			break;
		
		prev_v = v;				// Store the old vector and eigenvalue for the next iteration.
		prev_eigenvalue = eigenvalue;
	}
	
	if(v.getDotProduct(prev_v) < 0.0)			// Multiply eigenvalue with -1 if the dot product of the current vector and the previous vector is negative.
		eigenvalue *= -1.0;

	Matrix eigenvector = w / eigenvalue;
	
	return std::make_pair(eigenvalue, eigenvector);
}




// *** Main function ***
int main(int argc, char *argv[]) {
	try {
		int N = 0;								// N : Sizes of the A matrix, as A is a NxN square matrix.
		std::ifstream input_file(argv[1]);		// Input file stream declaration. (e.g. argv[1] = "A.txt")
		if(!input_file.is_open())				// Throw an exception if the input files with given filename does not exist.
			throw argv[1];
		
		// Count the number of rows in the A matrix with a temporary local string variable "line".
		for(std::string line; std::getline(input_file, line); N++);
		
		double **mat = new double*[N];			// Since we know N, now we can add the matrix declaration.
		
		input_file.clear();						// Clears any fail bits from the input file stream such as the end-of-file bit.
		input_file.seekg(0, input_file.beg);	// Seeks to the beginning of the input file.
		
		for(int i = 0; i < N; i++) {
			mat[i] = new double[N];				// Allocate memory for N integers to the i'th row of A matrix.
			for(int j = 0; j < N; j++)
				input_file >> mat[i][j];		// Get each entry of the A matrix from the input file stream.
		}
		
		input_file.close();						// Close the input file.
		
		Matrix A(N, N, mat);					// Create a NxN Matrix class object.
		
		double tol = atof(argv[2]);				// Since argv[2] is a char array (string), convert it to type double.
		
		// A.printMatrix("A");
		std::pair<double, Matrix> eigenval_and_eigenvec_1 = A.normPowIteration(tol);
		double eigenvalue_1 = eigenval_and_eigenvec_1.first;
		Matrix eigenvector_1 = eigenval_and_eigenvec_1.second;
		
		std::ofstream output_file(argv[3]);		// Output file stream declaration.
		if(!output_file.is_open())				// Throw an exception if the output file with given filename cannot be generated.
			throw argv[3];
		
		output_file << "Eigenvalue#1: " << eigenvalue_1 << std::endl;	// Print the dominant eigenvalue to the output file.
		
		for(int i = 0; i < N; i++)										// Print the corresponding dominant eigenvector to the output file.
			output_file << eigenvector_1.getEntry(i) << std::endl;
		
		Matrix deflated_A = A.getDeflatedMatrix(eigenvector_1, tol);	// Deflate the A matrix then calculate the second eigenvalue.
		double eigenvalue_2 = deflated_A.normPowIteration(tol).first;
		
		output_file << "Eigenvalue#2: " << eigenvalue_2;				// Print the second eigenvalue to the output file.
		
		output_file.close();					// Close the output file.
		
		return 0;				// Successful termination.
	}
	catch(const char *filename) {
		std::cerr << "[!ERROR!]: Could not open file \"" << filename << "\"." << std::endl;
		return 1;				// Unsuccessful termination.
	}
}

