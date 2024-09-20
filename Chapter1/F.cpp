#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>
using namespace std;

const double Pi = acos(-1);

class F1 : public Function{
public:
	double operator() (double x) const{
		return 89*sin(11.5/180*Pi)*sin(x)*cos(x)+89*cos(11.5/180*Pi)*pow(sin(x),2)-((49+55/2)*sin(11.5/180*Pi)-0.5*55*tan(11.5/180*Pi))*cos(x)-((49+55/2)*cos(11.5/180*Pi)-55/2)*sin(x);
	}
};

class F2 : public Function{
public:
	double operator() (double x) const{
		return 89*sin(11.5/180*Pi)*sin(x)*cos(x)+89*cos(11.5/180*Pi)*pow(sin(x),2)-((49+15)*sin(11.5/180*Pi)-15*tan(11.5/180*Pi))*cos(x)-((49+15)*cos(11.5/180*Pi)-15)*sin(x);
	}
};

void solve_a(){
	Newton_Method solver(F1(), Pi/4);
	double x = solver.solve();
	if(x==-1){cout<<"could not find root"<< endl;return;}
	cout<<"a root is: "<< x*180/Pi << endl;
	return;
}

void solve_b(){
	Newton_Method solver(F2(), 33*Pi/180);
	double x = solver.solve();
	if(x==-1){cout<<"could not find root"<< endl;return;}
	cout<<"a root is: "<< x*180/Pi << endl;
	return;
}

void solve_c(){
	Secant_Method solver_1(F2(), 33*Pi/180, 3*Pi);
	double x = solver_1.solve();
	cout<<"a root is: "<< x*180/Pi << endl;
	Secant_Method solver_2(F2(), 33*Pi/180, 2*Pi);
	x = solver_2.solve();
	cout<<"a root is: "<< x*180/Pi << endl;
	Secant_Method solver_3(F2(), 33*Pi/180, Pi);
	x = solver_3.solve();
	cout<<"a root is: "<< x*180/Pi << endl;
	return;
}

int main(){
	solve_a();
	solve_b();
	solve_c();
	return 0;
}
