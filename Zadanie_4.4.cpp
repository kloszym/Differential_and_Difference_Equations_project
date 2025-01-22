#include <iostream>
#include <vector>
#include <format>
#include <sstream>
#include "gnuplot-iostream.h"

typedef double (*EIFunction)(int, double, int);


int n;
double** B;
double* L;
EIFunction* eFunctions;
double epsilons[] = { sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)),sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)),-sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt(6.0 / 5.0)),-sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt(6.0 / 5.0)) };
double weights[] = { (18.0 - sqrt(30)) / 36.0,(18.0 + sqrt(30)) / 36.0,(18.0 + sqrt(30)) / 36.0,(18.0 - sqrt(30)) / 36.0 };
const double G = 0.000000000066743015;
const double roPi4G = 83.85;

std::string to_precise_string(double value) {
	std::ostringstream oss;
	oss << std::setprecision(15) << value;
	return oss.str();
}

void printing_matrix(double** matrix, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << matrix[i][j] << " ";
		}
		std::cout << "\n";
	}
}

double** allocMatrix2D(int n) {
	double** matrix2d = new double* [n];
	double* dumm = new double[n * n];
	for (int i = 0; i < n; i++) matrix2d[i] = dumm + i * n;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			matrix2d[i][j] = 0;
		}
	}
	return matrix2d;
}

double eI(int i, double x, int N) {
	double xI = 3.0 * i / N;
	double xIminus1 = 3.0 * (i - 1.0) / N;
	double xIplus1 = 3.0 * (i + 1.0) / N;
	if (x > xIminus1 && x <= xI) {
		return (x - xIminus1) / (xI - xIminus1);
	}
	else if (x > xI && x < xIplus1) {
		return (xIplus1 - x) / (xIplus1 - xI);
	}
	else {
		return 0.0;
	}
}
double eIPrim(int i, double x, int N) {
	double xI = 3.0 * i / N;
	double xIminus1 = 3.0 * (i - 1.0) / N;
	double xIplus1 = 3.0 * (i + 1.0) / N;
	if (x > xIminus1 && x <= xI) {
		return 1.0 / (xI - xIminus1);
	}
	else if (x > xI && x < xIplus1) {
		return -1.0 / (xIplus1 - xI);
	}
	else {
		return 0.0;
	}
}

double gaussianQuadratureOfOneFunction(EIFunction f, double a, double b, int j, int N) {
	// z powodu przybliżeń całki kwadraturą Gausaa możemy nie trafić w naszą funkcję
	// zatem lepiej wybrać granice całkowania tak by ograniczały się do przecięcia zbioru (a,b) i dziedziny funkcji e(j)
	double xJminus1 = 3.0 * (j - 1.0) / N;
	double xJplus1 = 3.0 * (j + 1.0) / N;

	if (xJminus1 > a) {

		a = xJminus1;

	}

	if (xJplus1 < b) {

		b = xJplus1;

	}

	double calka = 0;
	for (int i = 0; i < 4; i++) {
		calka = calka + ((f(j, ((a + b) / 2.0) + (((b - a) / 2.0) * epsilons[i]), N)) * weights[i]);
	}
	return calka * (b - a) / 2.0;
}

double gaussianQuadratureOfTwoFunctions(EIFunction f1, EIFunction f2, double a, double b, int j, int k, int N) {
	// Ten sam przypadek co w kwadraturze pojedynczej funkcji, jednakże a i b są brzegowmi wartościami naszej dziedziny, więc możemy to zapisać inaczej
	// całka oznaczona przyjmuje wartości różne od zera jedynie w
	// Zakładam, że j<=k
	double xJminus1 = 3.0 * (j - 1.0) / N;
	double xJplus1 = 3.0 * (j + 1.0) / N;
	double xKminus1 = 3.0 * (k - 1.0) / N;
	double xKplus1 = 3.0 * (k + 1.0) / N;

	a = xKminus1;
	b = xJplus1;

	double calka = 0;
	for (int i = 0; i < 4; i++) {
		calka = calka + ((f1(j, ((a + b) / 2.0) + (((b - a) / 2.0) * epsilons[i]), N) * f2(k, ((a + b) / 2.0) + (((b - a) / 2.0) * epsilons[i]), N)) * weights[i]);
	}
	return calka * (b - a) / 2.0;
}

double* solveEquation(double** b, double* l, int size) {
	int n = size;

	for (int p = 0; p < n; p++) {
		int max = p;
		for (int i = p + 1; i < n; i++) {
			if (abs(b[i][p]) > abs(b[max][p])) {
				max = i;
			}
		}
		double* temp = b[p];
		b[p] = b[max];
		b[max] = temp;

		double t = l[p];
		l[p] = l[max];
		l[max] = t;

		if (abs(b[p][p]) <= 0.00000000001) {
			std::cout << "blad";
			break;
		}

		for (int i = p + 1; i < n; i++) {
			double alpha = b[i][p] / b[p][p];
			l[i] -= alpha * l[p];
			for (int j = p; j < n; j++) {
				b[i][j] -= alpha * b[p][j];
			}
		}
	}


	double* x = new double[n];
	for (int i = n - 1; i >= 0; i--) {
		double sum = 0.0;
		for (int j = i + 1; j < n; j++) {
			sum += b[i][j] * x[j];
		}
		x[i] = (l[i] - sum) / b[i][i];
	}
	return x;
}



int main()
{
	Gnuplot gp;
	std::cout << "Na jak duzo przedzialow ma zostac podzielona dziedzina?\n";
	std::cin >> n;
	int ilePrzedzialow = n;
	n = n - 1;
	B = allocMatrix2D(n);

	L = new double[n];


	for (int i = 0; i < n; i++) {
		B[i][i] = gaussianQuadratureOfTwoFunctions(eIPrim, eIPrim, 0.0, 3.0, i + 1, i + 1, ilePrzedzialow);
	}
	for (int i = 0; i < n - 1; i++) {
		double answer = gaussianQuadratureOfTwoFunctions(eIPrim, eIPrim, 0.0, 3.0, i + 1, i + 2, ilePrzedzialow);
		B[i][i + 1] = answer;
		B[i + 1][i] = answer;
	}
	for (int i = 0; i < n; i++) {
		L[i] = gaussianQuadratureOfOneFunction(eI, 1, 2, i + 1, ilePrzedzialow) * roPi4G;
	}

	double* answOneToNMinus1 = solveEquation(B, L, n);
	double* answ = new double[n + 2];
	answ[0] = -5.0 * G;
	for (int i = 1; i < n + 2; i++) {
		answ[i] = answOneToNMinus1[i - 1];
	}
	answ[n + 1] = -4.0 * G;

	std::string przedzialy = std::to_string(ilePrzedzialow);

	std::string plotting = "plot $eI(0, x, " + przedzialy + ", " + to_precise_string(-5.0 * G) + ")";
	for (int i = 0; i < n; i++) {
		plotting = plotting + "+$eI(" + std::to_string(i + 1) + ", x, " + przedzialy + ", " + to_precise_string(answ[i]) + ")";
	}
	plotting = plotting + "+$eI(" + przedzialy + ", x, " + przedzialy + ", " + to_precise_string(-4.0 * G) + ") notitle\n";

	gp << "set terminal wxt\n";
	gp << "function $eI(i, x, N, mnoznik) << EOD\n";
	gp << "xI = 3.0 * i / N\n";
	gp << "xIminus1 = 3.0 * (i - 1.0) / N\n";
	gp << "xIplus1 = 3.0 * (i + 1.0) / N\n";
	gp << "if (x > xIminus1 && x <= xI) {\n";
	gp << "return mnoznik * (x - xIminus1) / (xI - xIminus1)}\n";
	gp << "else if (x > xI && x < xIplus1) {\n";
	gp << "return mnoznik * (xIplus1 - x) / (xIplus1 - xI)}\n";
	gp << "else {\n";
	gp << "return 0.0}\n";
	gp << "EOD\n";
	gp << "set xrange[0:3]\n";
	gp << plotting;
	gp << "pause -1\n";


	delete[] B[0];
	delete[] B;
	delete[] L;
	return 0;
}
