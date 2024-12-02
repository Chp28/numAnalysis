#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "spline.h"
using namespace std;

double F1(double x){
	return 1.0/(1+pow(x,2));
}

void test_A(){
	cout << "Problem A:"<< endl;
	vector<double> N = {6,11,21,41,81};
	vector<double> error_list;
	for(int i=0; i<5; i++){
		vector<double> x;
		vector<double> y;
		for(double j=-1; j<=1;){
			x.emplace_back(j);
			y.emplace_back(F1(j));
			j += 2/(N[i]-1);
		}
		double left = (F1(-1+1e-6)-F1(-1))/1e-6;
		double right = (F1(1)-F1(1-1e-6))/1e-6;
		PP_3_2_complete s(x, y, left, right);
		double max = 0;
		for(double j=-1+2/(N[i]-1); j<=1;){
			double error = abs(s.solve((j-1/(N[i]-1)))-F1((j-1/(N[i]-1))));
			if(max<error){max = error;}
			j += 2/(N[i]-1);
		}
		error_list.emplace_back(max);
	}
	for(int i=0; i<error_list.size(); i++){
		cout << error_list[i] << " ";
	}
	cout << endl;
	return;
}

void test_D(){
	cout << "Problem D:" << endl;
	vector<double> x;
	vector<double> y;
	for(double i=-5; i<=5; i++){
		x.emplace_back(i);
		y.emplace_back(F1(i));
	}
	Interpolation_BSpline s(x,y,3);
	s.complete_Spline_3_2((F1(-5+1e-6)-F1(-5))/1e-6,(F1(5)-F1(5-1e-6))/1e-6);
	vector<double> temp = {-3.5,-3,-0.5,0,0.5,3,3.5};
	for(int i=0; i<=6; i++){
		cout << abs(s.solve(temp[i])-F1(temp[i])) << " ";
	}
	cout << endl;
	MatrixXf A = MatrixXf::Zero(10,10);
	VectorXf B = VectorXf::Zero(10);
	VectorXf t = VectorXf::Zero(10);
	for(int i=0; i<10; i++){
		if(i==0){
			A(i,i) = 5;
			A(i,i+1) = 1;
			B(i) = 8*F1(-4.5)-2*F1(-5);
		}else{
			if(i==9){
				A(i,i-1) = 1;
				A(i,i) = 5;
				B(i) = 8*F1(4.5)-2*F1(5);
			}else{
				A(i,i) = 6;
				A(i,i-1) = 1;
				A(i,i+1) = 1;
				B(i) = 8*F1(-4.5+i);
			}
		}
	}
	t = A.colPivHouseholderQr().solve(B);
	vector<double> coef;
	coef.emplace_back(2*F1(-5)-t[0]);
	for(int i=0; i<10; i++){
		coef.emplace_back(t[i]);
	}
	coef.emplace_back(2*F1(5)-t[9]);
	Interpolation_BSpline _s(x,2,coef);
	for(int i=0; i<=6; i++){
		cout << abs(_s.solve(temp[i])-F1(temp[i])) << " ";
	}
	cout << endl;
	
}

int main(){
	test_A();
	test_D();
	return 0;
}


