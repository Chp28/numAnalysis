#ifndef EQUATIONSOLVER
#define EQUATIONSOLVER

#include "Function.hpp"
#include <cmath>
#include <iostream>
using namespace std;

class EquationSolver{
protected:
	const Function & F;
public:
	EquationSolver(const Function& F) : F(F) {}
	virtual double solve() = 0;
};

class Bisection_Method : public EquationSolver {
private:
	double a,b;
	double eps, delta;
	int Maxiter;
public:
	Bisection_Method(const Function &F,double a,double b,double eps=1e-7,double delta=1e-6,int Maxiter=50) : EquationSolver(F), a(a), b(b), eps(eps), delta(delta), Maxiter(Maxiter) {}
	
	double solve() override{
		int k = 0;
		double c;
		for(k = 0; k < Maxiter; k++){
			if(F(a)*F(b)>0){
				return -1;
			}
			c = a+(b-a)/2;
			if((b-a)<delta){break;}
			if(abs(F(c))<eps){break;}
			else if(F(c)*F(a)<0){
				b = c;
			}else{
				a = c;
			}
		}
		if(k==Maxiter){
			return -1;
		}
		return c;
	}
};

class Newton_Method : public EquationSolver{
private:
	double x0;
	double eps;
	int Maxiter;
public:
	Newton_Method(const Function &F,double x0,double eps=1e-7,double Maxiter=8) :
	EquationSolver(F), x0(x0), Maxiter(Maxiter), eps(eps) {}
	
	double solve() override{
		int k = 0;
		for(k; k<Maxiter ; k++){
			if(F.derivative(x0)==0){
				return -1;
			}
			x0 = x0 - F(x0)/F.derivative(x0);
			if(abs(F(x0))<eps){
				return x0;
			}
		}
		return -1;
	}
};

class Secant_Method : public EquationSolver{
private:
	double x0,x1;
	double eps,delta;
	int Maxiter;
public:
	Secant_Method(const Function &F,double x0,double x1,double eps=1e-7,double delta=1e-6, double Maxiter=50) : EquationSolver(F), x0(x0), x1(x1), eps(eps), delta(delta), Maxiter(Maxiter) {}
	double solve() override{
		int k = 0;
		double s;
		for(k; k<Maxiter ; k++){
			s = x1;
			x1 = x1-F(x1)*(x1-x0)/(F(x1)-F(x0));
			x0 = s;
			if(abs(x1-x0)<delta){
				return x1;
			}else if(abs(F(x1))<eps){
				return x1;
			}
		}
		return -1;
	}
};

#endif

