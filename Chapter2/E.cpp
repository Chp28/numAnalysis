#include "Function.hpp"
#include "Interpolation.hpp"
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

double compute_error(vector<double> x,vector<double> y,Newton_formula solver){
	double sum = 0;
	for(int i=0;i<x.size();i++){
		sum += fabs(y[i]-solver.solve(x[i]));
	}
	return sum;
}

int main(){
	vector<double> x = {0,6,10,13,17,20,28};
	vector<double> y_1 = {6.67,17.3,42.7,37.3,30.1,29.3,28.7};
	vector<double> y_2 = {6.67,16.1,18.9,15.0,10.6,9.44,8.89};
	vector<double> best_1;
	vector<double> best_2;
	vector<double> temp;
	double Min_1 = 10000;
	double Min_2 = 10000;
	for(int i=0;i<7;i++){
		temp.emplace_back(i);
		for(int j=i+1;j<7;j++){
			temp.emplace_back(j);
			for(int k=j+1;k<7;k++){
				temp.emplace_back(k);
				vector<double> _x = {x[i] , x[j] , x[k]};
				vector<double> _y_1 = {y_1[i] , y_1[j] , y_1[k]};
				vector<double> _y_2 = {y_2[i] , y_2[j] , y_2[k]};
				Newton_formula solver_1(_x,_y_1);
				Newton_formula solver_2(_x,_y_2);
				solver_1.Generate();
				solver_2.Generate();
				double sum_1 = compute_error(x,y_1,solver_1);
				double sum_2 = compute_error(x,y_2,solver_2);
				if(Min_1>sum_1){
					Min_1 = sum_1;
					best_1 = temp;
				}
				if(Min_2>sum_2){
					Min_2 = sum_2;
					best_2 = temp;
				}
				temp.pop_back();
			}
			temp.pop_back();
		}
		temp.pop_back();
	}
	
	vector<double> _x_1 = {x[best_1[0]] , x[best_1[1]] , x[best_1[2]]};
	vector<double> _y_1 = {y_1[best_1[0]] , y_1[best_1[1]] , y_1[best_1[2]]};
	vector<double> _x_2 = {x[best_2[0]] , x[best_2[1]] , x[best_2[2]]};
	vector<double> _y_2 = {y_2[best_2[0]] , y_2[best_2[1]] , y_2[best_2[2]]};
	Newton_formula solver_1(_x_1,_y_1);
	Newton_formula solver_2(_x_2,_y_2);
	solver_1.Generate();
	solver_2.Generate();
	
	cout << "Sp1 will choose " << best_1[0] << " " << best_1[1] << " " << best_1[2] << " to build the newton_formula" << endl;
	cout << "Sp2 will choose " << best_2[0] << " " << best_2[1] << " " << best_2[2] << " to build the newton_formula" << endl;
	
	solver_1.print_coef();
	solver_2.print_coef();
	
	if(solver_1.solve(43)<=0){
		cout << "Sp1 will die after another 15 days" << endl;
	}else{
		cout << "Sp1 will not die after another 15 days" << endl;
	}
	if(solver_2.solve(43)<=0){
		cout << "Sp2 will die after another 15 days" << endl;
	}else{
		cout << "Sp2 will not die after another 15 days" << endl;
	}
	
	return 0;
}
