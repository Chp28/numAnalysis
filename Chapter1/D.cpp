#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>
#define e 2.718281828459
using namespace std;

const double Pi = acos(-1);

class F1 : public Function{
public:
	double operator() (double x) const {
		return sin(x/2)-1;
	}
};

class F2 : public Function{
public:
	double operator() (double x) const {
		return pow(e,x)-tan(x);
	}
};

class F3 : public Function{
public:
	double operator() (double x) const {
		return pow(x,3)-12*pow(x,2)+3*x+1;
	}
};

void solve_f1(){
	cout << "Solving sin(x/2)-1 with x0 = 0, x1 = Pi/2" <<endl;
	Secant_Method solver_f1(F1(), 0, Pi/2);
	double x = solver_f1.solve();
	if(x == -1){
		cout << "could not find root" << endl;
		return;
	}
	cout << "A root is: " << x << endl;
	return;
}

void solve_f2(){
	cout << "Solving e^{x}-tan(x) with x0 = 1, x1 = 1.4" <<endl;
	Secant_Method solver_f2(F2(), 1, 1.4);
	double x = solver_f2.solve();
	if(x == -1){
		cout << "could not find root" << endl;
		return;
	}
	cout << "A root is: " << x << endl;
	return;
}

void solve_f3(){
	cout << "Solving x^{3}-12x^{2}+3x+1 with x0 = 0, x1 = -0.5" <<endl;
	Secant_Method solver_f3(F3(), 0, -0.5);
	double x = solver_f3.solve();
	if(x == -1){
		cout << "could not find root" << endl;
		return;
	}
	cout << "A root is: " << x << endl;
	return;
}

int main(){
	solve_f1();
	solve_f2();
	solve_f3();
	return 0;
}


