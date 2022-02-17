// *** Header Imports ***
#include <iostream>			// Import cout.
//#include <iomanip>		// Import setf() and setprecision() functions.
#include <cmath>			// Import fabs() and pow() function.
#include <utility>			// Import pair.

#define ITER_MAX 10000				// Maximum number of iterations allowed in each method.
#define BISECTION_ITER_MAX 2		// Maximum number of iterations allowed for the Bisection Method in the Hybrid Method.

// *** Function Prototypes ***
std::pair<int, double> bisectionMethod(double *coeff_arr, const int size, double x_0, double x_1, const double tol);
std::pair<int, double> secantMethod(double *coeff_arr, const int size, double x_0, double x_1, const double tol);
std::pair<int, double> hybridMethod(double *coeff_arr, const int size, double x_0, double x_1, const double tol);
double getYvalue(double *coeff_arr, const int size, const double x);
template<typename T> int sgn(const T val);

// *** Main function ***

// argc : Number of command line arguments, argv : Command line arguments as character arrays.
int main(int argc, char *argv[]) {
	// Command line arguments are: test a_n a_(n-1) ... a_0 x_0 x_1 tol
	// Since n >= 1, the number of command line arguments has to be at least 5 for the program to continue.
	if(argc < 5) {
		std::cerr << "[!ERROR!]: Number of command line arguments \"argc = " << argc << "\" has to be at least 5." << std::endl;
		return 1;				// Unsuccessful termination.
	}
	
	double *coeff_arr = new double[argc - 4];	// Create a dynamically allocated coefficient array for the polynomial where coeff_arr[i] is equal to a_(n-i)
	for(int i = 1; i < argc - 3; i++)			// argv[0] is the executable filename so we store the command line arguments 1, 2, ..., argc - 4
		coeff_arr[i-1] = atof(argv[i]);
	
	double x_0 = atof(argv[argc - 3]);			// Remaining three command line arguments are x_0, x_1 and tolerance respectively.
	double x_1 = atof(argv[argc - 2]);			// We use the atof function to convert character arrays (strings) to floating point numbers.
	double tol = atof(argv[argc - 1]);
	
	if(x_0 >= x_1) {							// x_0 < x_1 has to hold for the program to continue.
		std::cerr << "[!ERROR!]: x_0 = " << x_0 << "  >=  x_1 = " << x_1 << std::endl;
		return 1;				// Unsuccessful termination.
	}
	
	if(tol <= 0.0) {							// Tolerance value has to be positive for the program to continue.
		std::cerr << "[!ERROR!]: Tolerance value tol = " << tol << " has to be positive. " << std::endl;
		return 1;				// Unsuccessful termination.
	}
	
	// First number is the number of iterations, second number is the root for each pair.
	
	// Print the calculated root and the number of iterations for the Bisection Method.
	std::pair<int, double> sol1 = bisectionMethod(coeff_arr, argc - 4, x_0, x_1, tol);
	std::cout << "For the Bisection Method: " << sol1.first << " number of iterations and " << sol1.second << " is the root." << std::endl;
	
	// Print the calculated root and the number of iterations for the Secant Method.
	std::pair<int, double> sol2 = secantMethod(coeff_arr, argc - 4, x_0, x_1, tol);
	std::cout << "For the Secant Method: " << sol2.first << " number of iterations and " << sol2.second << " is the root." << std::endl;
	
	// Print the calculated root and the number of iterations for the Hybrid Method.
	std::pair<int, double> sol3 = hybridMethod(coeff_arr, argc - 4, x_0, x_1, tol);
	std::cout << "For the Hybrid Method: " << sol3.first << " number of iterations and " << sol3.second << " is the root." << std::endl;
	
	delete [] coeff_arr;	// Deallocate the coefficient array.
	return 0;				// Successful termination.
}


// *** Function Implementations ***

/**
 * Implements the bisection method to calculate the root of the polynomial within the given tolerance value.
 * @param coeff_arr Coefficient array of the polynomial where coefficients are ordered from a_n to a_0.
 * @param size Size of the coefficient array, which is basically n+1 if the polynomial is of order n.
 * @param x_0 Lower initial guess.
 * @param x_1 Upper initial guess.
 * @param tol Tolerance value.
 * @return Number of iterations and the calculated root of the polynomial respectively.
 *
 * Pseudocode from the book:
 *
 * while ((b-a) > tol) do
 *     m = a + (b - a) / 2
 *     if sign(f(a)) = sign(f(m)) then
 *         a = m
 *     else
 *         b = m
 *     end
 * end
 */
std::pair<int, double> bisectionMethod(double *coeff_arr, const int size, double x_0, double x_1, const double tol) {
	int iteration_count = 1;				// A counter which keeps track of the iteration number, starts from 1 and increments by 1 at the end of the loop.
	double root;
	
	while(iteration_count < ITER_MAX) {		// If the iteration limit is reached then there is no point in doing anymore iterations, so break out of the loop.
		root = x_0 + (x_1 - x_0) / 2.0;		// (a + b) / 2 would work too but a + (b - a) / 2 prevents overflow.
		
		if(sgn(getYvalue(coeff_arr, size, x_0)) == sgn(getYvalue(coeff_arr, size, root)))
			x_0 = root;		// If f(a) and f(m) both have the same signs, we can choose the next lower bound to be a = m.
		else
			x_1 = root;		// If their signs differ, we choose the next upper bound as b = m.
		
		if(std::fabs(x_1 - x_0) <= tol)		// If abs(b - a) is lower than the tolerance value then terminate the process.
			break;
		
		iteration_count++;
	}
	
	return std::make_pair(iteration_count, root);
}

/**
 * Implements the secant method to calculate the root of the polynomial within the given tolerance value.
 * @param coeff_arr Coefficient array of the polynomial where coefficients are ordered from a_n to a_0.
 * @param size Size of the coefficient array, which is basically n+1 if the polynomial is of order n.
 * @param x_0 Lower initial guess.
 * @param x_1 Upper initial guess.
 * @param tol Tolerance value.
 * @return Number of iterations and the calculated root of the polynomial respectively.
 *
 * Pseudocode from the book:
 *
 * x_0, x_1 = initial guesses
 * for k = 1, 2, ...
 *     x_(k+1) = x_k - f(x_k) * (x_k - x_(k-1)) / (f(x_k) - f(x_(k-1)))
 * end
 */
std::pair<int, double> secantMethod(double *coeff_arr, const int size, double x_0, double x_1, const double tol) {
	int iteration_count = 1;				// A counter which keeps track of the iteration number, starts from 1 and increments by 1 at the end of the loop.
	double root;
	
	while(iteration_count < ITER_MAX) {		// If the iteration limit is reached then there is no point in doing anymore iterations, so break out of the loop.
		const double val_0 = getYvalue(coeff_arr, size, x_0);	// val_0 : f(x_(k-1))
		const double val_1 = getYvalue(coeff_arr, size, x_1);	// val_1 : f(x_k)
		root = x_1 - val_1 * (x_1 - x_0) / (val_1 - val_0);		// root  : x_(k+1)
		
		if(std::fabs(root - x_1) <= tol)	// If abs(x_(k+1) - x_k) is lower than the tolerance value then terminate the process.
			break;
		
		x_0 = x_1;				// Shift the x values by one for the next iteration. This way we preserve memory.
		x_1 = root;				// Before :    x_0 : x_(k-1), x_1 : x_k,     root : x_(k+1)
		iteration_count++;		// After  :    x_0 : x_k,     x_1 : x_(k+1), root : x_(k+1)
	}
	
	return std::make_pair(iteration_count, root);
}

/**
 * Combines both the bisection and the secant method to calculate the root of the polynomial within the given tolerance value.
 * @param coeff_arr Coefficient array of the polynomial where coefficients are ordered from a_n to a_0.
 * @param size Size of the coefficient array, which is basically n+1 if the polynomial is of order n.
 * @param x_0 Lower initial guess.
 * @param x_1 Upper initial guess.
 * @param tol Tolerance value.
 * @return Number of iterations and the calculated root of the polynomial respectively.
 *
 * Pseudocode:
 *
 * k = 1
 *
 * while ((b-a) > tol and k <= 2) do
 *     m = a + (b - a) / 2
 *     if sign(f(a)) = sign(f(m)) then
 *         a = m
 *     else
 *         b = m
 *     end
 *     k = k + 1
 * end
 *
 * x_2 = a
 * x_3 = b
 * for k = 3, 4, ...
 *     x_(k+1) = x_k - f(x_k) * (x_k - x_(k-1)) / (f(x_k) - f(x_(k-1)))
 * end
 */
std::pair<int, double> hybridMethod(double *coeff_arr, const int size, double x_0, double x_1, const double tol) {
	int iteration_count = 1;				// A counter which keeps track of the iteration number, starts from 1 and increments by 1 at the end of the loop.
	double root;
	
	while(iteration_count <= BISECTION_ITER_MAX) {			// The Bisection Method is used for the first 2 iterations.
		root = x_0 + (x_1 - x_0) / 2.0;		// (a + b) / 2 would work too but a + (b - a) / 2 prevents overflow.
		
		if(sgn(getYvalue(coeff_arr, size, x_0)) == sgn(getYvalue(coeff_arr, size, root)))
			x_0 = root;		// If f(a) and f(m) both have the same signs, we can choose the next lower bound to be a = m.
		else
			x_1 = root;		// If their signs differ, we choose the next upper bound as b = m.
		
		if(std::fabs(x_1 - x_0) <= tol)		// If abs(x_(k+1) - x_k) is lower than the tolerance value then terminate the process and return the calculated root along with the iteration number.
			return std::make_pair(iteration_count, root);
		
		iteration_count++;
	}
	
	while(iteration_count < ITER_MAX) {						// For the rest of the iterations the Secant Method is used.
		const double val_0 = getYvalue(coeff_arr, size, x_0);	// val_0 : f(x_(k-1))
		const double val_1 = getYvalue(coeff_arr, size, x_1);	// val_1 : f(x_k)
		root = x_1 - val_1 * (x_1 - x_0) / (val_1 - val_0);		// root  : x_(k+1)
		
		if(std::fabs(root - x_1) <= tol)	// If abs(x_(k+1) - x_k) is lower than the tolerance value then break and terminate the process.
			break;
		
		x_0 = x_1;				// Shift the x values by one for the next iteration. This way we preserve memory.
		x_1 = root;				// Before :    x_0 : x_(k-1), x_1 : x_k,     root : x_(k+1)
		iteration_count++;		// After  :    x_0 : x_k,     x_1 : x_(k+1), root : x_(k+1)
	}
	
	return std::make_pair(iteration_count, root);
}

/**
 * Calculates and returns the y = f(x) value for the given x.
 * @param coeff_arr Coefficient array of the polynomial where coefficients are ordered from a_n to a_0.
 * @param size Size of the coefficient array, which is basically n+1 if the polynomial is of order n.
 * @param x The given x value.
 * @return f(x) for the given x value.
 */
double getYvalue(double *coeff_arr, const int size, const double x) {
	double y = 0.0;
	for(int i = 0; i < size; i++)
		y += coeff_arr[i] * std::pow(x, size - i - 1);		// Add each term.
	return y;	// y = a_n * x^n + a_(n-1) * x^(n-1) + ... + a_1 * x + a_0
}

/**
 * A type-safe template function which returns the sign of the given number.
 * @param val The given number.
 * @return 1 if val is positive, -1 if val is negative, 0 otherwise.
 */
template<typename T> int sgn(const T val) {
	return (T(0) < val) - (val < T(0));
}
