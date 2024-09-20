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
		return 1.0/x-tan(x);
	}
};

class F2 : public Function{
public:
	double operator() (double x) const{
		return 1.0/x-pow(2,x);
	}
};

class F3 : public Function{
public:
	double operator() (double x) const{
		return pow(2,-x)+pow(e,x)+2*cos(x)-6;
	}
};

class F4 : public Function{
public:
	double operator() (double x) const{
		return (pow(x,3)+4*pow(x,2)+3*x+5)/(2*pow(x,3)-9*pow(x,2)+18*x-2);
	}
};

void solve_f1(){
	cout << "Solving x^{-1} - tan(x) on [0,Pi/2]" <<endl;
	Bisection_Method solver_f1(F1(), 0, Pi/2);
	double x = solver_f1.solve();
	if(x == -1){
		cout << "could not find root" << endl;
		return;
	}
	cout << "A root is: " << x << endl;
	return;
}

void solve_f2(){
	cout << "Solving x^{-1}-2^{x} on [0,1]" <<endl;
	Bisection_Method solver_f2(F2(), 0, 1);
	double x = solver_f2.solve();
	if(x == -1){
		cout << "could not find root" << endl;
		return;
	}
	cout << "A root is: " << x << endl;
	return;
}

void solve_f3(){
	cout << "Solving 2^{-x}+e^{x}+2cos(x)-6 on [1,3]" <<endl;
	Bisection_Method solver_f3(F3(), 1, 3);
	double x = solver_f3.solve();
	if(x == -1){
		cout << "could not find root" << endl;
		return;
	}
	cout << "A root is: " << x << endl;
	return;
}

void solve_f4(){
	cout << "Solving (x^{3}+4x^{2}+3x+5)/(2x^{3}-9x^{2}+18x-2) on [0,4]" <<endl;
	Bisection_Method solver_f4(F4(), 0, 4);
	double x = solver_f4.solve();
	if(x == -1){
		cout << "no root" << endl;
		return;
	}
	cout << "A root is: " << x << endl;
	return;
}


int main(){
	solve_f1();
	solve_f2();
	solve_f3();
	solve_f4();
	return 0;
}
