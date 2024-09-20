#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>
using namespace std;

const double Pi = acos(-1);

class F : public Function{
public:
	double operator() (double x) const {
		return 10*(0.5*Pi-asin(x)-x*pow(1-pow(x,2),0.5))-12.4;
	}
};

void solve_f1(){
	cout << "use bisection method to solve 10(Pi/2-asin(x)-x(1-x^{2})^{1/2})" <<endl;
	Bisection_Method solver_f1(F(), 0, 1);
	double x = solver_f1.solve();
	if(x == -1){
		cout << "could not find root" << endl;
		return;
	}
	cout << "A root is: " << x << endl;
	return;
}

void solve_f2(){
	cout << "use Newton method to solve 10(Pi/2-asin(x)-x(1-x^{2})^{1/2})" <<endl;
	Newton_Method solver_f2(F(), 0.1);
	double x = solver_f2.solve();
	if(x == -1){
		cout << "could not find root" << endl;
		return;
	}
	cout << "A root is: " << x << endl;
	return;
}

void solve_f3(){
	cout << "use bisection method to solve 10(Pi/2-asin(x)-x(1-x^{2})^{1/2})" <<endl;
	Secant_Method solver_f3(F(), 0, 1);
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
