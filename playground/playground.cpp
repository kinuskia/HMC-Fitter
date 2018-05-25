#include <iostream>
#include <cmath>
#include "../auxiliary_files/vector.hpp"
#include "../auxiliary_files/matrix.hpp"
#include <random>
#include <assert.h>

typedef double number_type;
typedef std::size_t size_type;






int main ()
{

	Matrix<number_type> A(3, 3);
	A[0][0] = 6.;
	A[0][1] = 2.;
	A[0][2] = 1.;
	A[1][0] = -3.;
	A[1][1] = 1.;
	A[1][2] = 6.;
	A[2][0] = -6.;
	A[2][1] = 8.;
	A[2][2] = 10.;
	Matrix<number_type> A_inv(A.rowsize(), A.colsize());
	MatrixInverse(A, A_inv);

	Vector<number_type> b(3);
	b[0] = -2.;
	b[1] = 1.;
	b[2] = 3.;

	Vector<number_type> x(3);

	Vector<number_type> s(A.colsize());
	Vector<size_type> p(A.colsize());
	Vector<size_type> q(A.colsize());

	row_equilibrate(A, s);
	lu_fullpivot(A, p, q);
	x = number_type(0.0);
	apply_equilibrate(s, b);
	permute_forward(p, b);
	solveL(A, b, b);
	solveU(A, x, b);
	permute_backward(q, x);

	for (size_type i = 0; i < x.size(); ++i)
		std::cout << x[i] << "\n";

	for (size_type i = 0; i < A_inv.rowsize(); ++i)
	{
		for (size_type j = 0; j < A_inv.colsize(); ++j)
		{
			std::cout << A_inv(i, j) << " ";
		}
		std::cout << "\n";
	}


	return 0;
}