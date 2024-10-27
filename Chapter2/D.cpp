#include "Function.hpp"
#include "Interpolation.hpp"
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

int main(){
	vector<double> x = {0,3,5,8,13};
	vector<double> y = {0,225,383,623,993};
	vector<double> x_d = {75,77,80,74,72};
	Hermite solver(x,y,x_d);
	cout << "when t = 10 , displacement is " << solver.solve(10) << endl;
	cout << "when t = 10 , velocity is " << solver.solve_derivative(10) << endl;
	
	for(double k=0; k<13; k+=0.1){
		double d = solver.solve_derivative(k);
		if(d > 81){
			cout << "velocity is " << d << " at " << k << endl;
		}
	}
	return 0;
}
