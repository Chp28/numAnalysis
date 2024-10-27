#ifndef INTERPOLATION
#define INTERPOLATION

#include "Function.hpp"
#include <vector> 
#include <iostream>
#include <cmath>
using namespace std;

class useless : public Function{
public:
	double operator()(double x) const{
		return 0;
	}
};

class Interpolation{
protected:
	const vector<double>& x;
public:
	Interpolation(const vector<double>& x) : x(x) {}
	virtual vector<double> solve_vector(const vector<double>& z) = 0;
	virtual double solve_derivative(const double& z) = 0;
	virtual double solve(const double& z) = 0;
};

class Newton_formula : public Interpolation {
protected:
	const Function& F;
	vector<double> y;
	vector<vector<double>> coef;
public:
	Newton_formula(const vector<double>& x, vector<double> y) : F(useless()) , Interpolation(x) , y(y) {}
	
	Newton_formula(const Function& F, const vector<double>& x) : F(F) , Interpolation(x)  {
		int size = x.size();
		for(int i = 0; i<size ; i++){
			y.emplace_back(F(x[i]));
		}
	}
	
	void Generate(){
		int size = x.size();
		coef.resize(size);
		int i,j;
		for(i=0; i<size; i++){
			for(j=0; j<size; j++){
				if(i==0){
					coef[i].emplace_back(y[j]);
				}else{
					coef[i].emplace_back(0);
				}
			}
		}
		for(i=1; i<size; i++){
			for(j=i; j<size; j++){
				coef[i][j] = (coef[i-1][j]-coef[i-1][j-1])/(x[j]-x[j-i]);
			}
		}
	}
	
	double solve(const double& z) override{
		double z_y = coef[0][0];
		double a = 1;
		for(int j=1; j<coef.size(); j++){
			a = a*(z-x[j-1]);
			z_y += a*coef[j][j];
		}
		return z_y;
	}
	
	vector<double> solve_vector(const vector<double>& z) override{
		int size = z.size();
		vector<double> z_y;
		for(int i=0; i<size; i++){
			z_y.emplace_back(solve(z[i]));
		}
		return z_y;
	}
	
	double solve_derivative(const double& z) override{
		return (solve(z+1e-6)-solve(z))/1e-6;
	}
	
	void print_coef(){
		if(coef.size()==0){
			cout << "have not yet generated coefficient" << endl;
		}else{
			cout << "The coefficient are ";
			for(int i=0; i<coef.size(); i++){
				cout << coef[i][i] << " ";
			}
			cout << "in order" << endl;
		}
	}
	
};

class Hermite : public Interpolation {
protected:
	const vector<double>& y;
	const vector<double>& x_d;
public:
	Hermite(const vector<double>& x,const vector<double>& y, const vector<double>& x_d) : Interpolation(x) , y(y) , x_d(x_d) {}
	
	double value(int i, const double& z){
		vector<double> k(4);
		int j;
		double sum;
		for(j=0; j<4; j++){
			sum = pow(1-(z-x[i])/(x[i+1]-x[i]),3-j)*pow((z-x[i])/(x[i+1]-x[i]),j);
			k[j] = sum;
		}
		return y[i]*k[0]+(y[i]+((double) 1/3)*x_d[i]*(x[i+1]-x[i]))*3*k[1]+(y[i+1]-((double) 1/3)*x_d[i+1]*(x[i+1]-x[i]))*3*k[2]+y[i+1]*k[3];
	}
	
	double solve(const double& z) override{
		int i,j;
		int size = x.size();
		double z_y;
		for(i=0; i<size-1 ; i++){
			if(z>=x[i] && z<=x[i+1]){
				z_y = value(i,z);
				break;
			}
		}
		return z_y;
	}
	
	vector<double> solve_vector(const vector<double>& z) override{
		int size = z.size();
		vector<double> z_y;
		for(int i=0; i<size ; i++){
			z_y.emplace_back(solve(z[i]));
		}
		return z_y;
	}
	
	double solve_derivative(const double& z){
		return (solve(z+1e-8)-solve(z))/1e-8;
	}
};

#endif
