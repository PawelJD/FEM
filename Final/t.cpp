#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>

#include "t.h"

using std::vector;
using std::cout;
using std::endl;

static void luDecompose(vector<vector<double>>& A, vector<int>& piv)
{
	int n = (int)A.size();
	piv.resize(n);
	for (int i = 0; i < n; ++i) piv[i] = i;

	for (int k = 0; k < n; ++k) {
		// pivot
		int maxRow = k;
		double maxVal = std::abs(A[k][k]);
		for (int i = k + 1; i < n; ++i) {
			double v = std::abs(A[i][k]);
			if (v > maxVal) { maxVal = v; maxRow = i; }
		}

		if (maxVal < 1e-14) {
			std::cerr << "Blad: zerowy pivot w kolumnie " << k << "\n";
		}

		if (maxRow != k) {
			std::swap(A[k], A[maxRow]);
			std::swap(piv[k], piv[maxRow]);
		}

		double Akk = A[k][k];
		if (std::abs(Akk) < 1e-14) continue;

		for (int i = k + 1; i < n; ++i) {
			double factor = A[i][k] / Akk;
			A[i][k] = factor;
			for (int j = k + 1; j < n; ++j) {
				A[i][j] -= factor * A[k][j];
			}
		}
	}
}

static void luSolve(const vector<vector<double>>& LU,
	const vector<int>& piv,
	const vector<double>& b,
	vector<double>& x)
{
	int n = (int)LU.size();
	x.assign(n, 0.0);

	// Pb
	vector<double> bp(n);
	for (int i = 0; i < n; ++i)
		bp[i] = b[piv[i]];

	vector<double> y(n, 0.0);
	for (int i = 0; i < n; ++i) {
		double sum = bp[i];
		for (int j = 0; j < i; ++j)
			sum -= LU[i][j] * y[j];
		y[i] = sum;
	}

	for (int i = n - 1; i >= 0; --i) {
		double sum = y[i];
		for (int j = i + 1; j < n; ++j)
			sum -= LU[i][j] * x[j];
		x[i] = sum / LU[i][i];
	}
}

void SolveT(Equation& eq, int n, double InitialTemp,
	int SimulationStepTime, int SimulationTime, int Tot)
{
	const double dt = (double)SimulationStepTime;
	const int steps = SimulationTime / SimulationStepTime;

	eq.resetT(n);

	vector<double> D(n, InitialTemp);

	vector<double> Tmin(steps, (double)Tot);
	vector<double> Tmax(steps, 0.0);

	vector<vector<double>> Cdt(n, vector<double>(n, 0.0));
	vector<vector<double>> A(n, vector<double>(n, 0.0));

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			Cdt[i][j] = eq.CG[i][j] / dt;
			A[i][j] = eq.HG[i][j] + Cdt[i][j];
		}
	}

	vector<int> piv;
	luDecompose(A, piv);

	vector<double> B(n, 0.0);
	vector<double> CT0(n, 0.0);
	vector<double> x; // rozwi¹zanie

	for (int l = 0; l < steps; ++l) {
		// CT0 = Cdt * D
		std::fill(CT0.begin(), CT0.end(), 0.0);
		for (int i = 0; i < n; ++i) {
			double s = 0.0;
			for (int j = 0; j < n; ++j)
				s += Cdt[i][j] * D[j];
			CT0[i] = s;
		}

		// B = PG + CT0
		for (int i = 0; i < n; ++i)
			B[i] = eq.PG[i] + CT0[i];

		// Rozwi¹¿ A * t = B
		luSolve(A, piv, B, x);

		for (int i = 0; i < n; ++i) {
			eq.t[i] = x[i];
			D[i] = x[i];

			if (x[i] < Tmin[l]) Tmin[l] = x[i];
			if (x[i] > Tmax[l]) Tmax[l] = x[i];
		}
	}

	cout << "\n\nMin i Max\n";
	for (int i = 0; i < steps; ++i) {
		cout << "Tmin = " << std::fixed << std::setprecision(6) << Tmin[i]
			<< ", Tmax = " << std::fixed << std::setprecision(6) << Tmax[i] << "\n";
	}
}


//void SolveT(Equation& eq, int n, double InitialTemp, int SimulationStepTime, int SimulationTime, int Tot) {
//	vector<vector<double>> A(n, vector<double>(n));
//	vector<double> B(n);
//
//	// Kopiowanie macierzy HG do A i wektora PG do B
//	for (int i = 0; i < n; ++i) {
//		for (int j = 0; j < n; ++j) {
//			A[i][j] = eq.HG[i][j];
//		}
//		B[i] = eq.PG[i];
//	}
//
//	// Eliminacja Gaussa
//	for (int k = 0; k < n; ++k) {
//
//		int maxRow = k;
//
//		for (int i = k + 1; i < n; ++i) {
//			if (abs(A[i][k]) > abs(A[maxRow][k])) {
//				maxRow = i;
//			}
//		}
//
//		if (fabs(A[maxRow][k]) < 1e-12) {
//			cerr << "Blad: zerowy pivot w kolumnie " << k << endl;
//			continue;
//		}
//
//		if (maxRow != k) {
//			for (int j = 0; j < n; j++) {
//				swap(A[k][j], A[maxRow][j]);
//			}
//			swap(B[k], B[maxRow]);
//		}
//
//		for (int i = k + 1; i < n; ++i) {
//			double factor = A[i][k] / A[k][k];
//			for (int j = k; j < n; ++j) {
//				A[i][j] -= factor * A[k][j];
//			}
//			B[i] -= factor * B[k];
//		}
//
//	}
//
//	eq.resetT(n);
//
//	for (int i = n - 1; i >= 0; --i) {
//		double sum = B[i];
//		for (int j = i + 1; j < n; ++j) {
//			sum -= A[i][j] * eq.t[j];
//		}
//		eq.t[i] = sum / A[i][i];
//	}
//
//
//}