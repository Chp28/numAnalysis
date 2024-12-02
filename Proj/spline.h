#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
using namespace std;
using namespace Eigen;

class B_spline{
public:
	B_spline(vector<double> point, int N) : N(N){
		knots.emplace_back(point[0]-1);
		for(int j=0; j<point.size(); j++){
			knots.emplace_back(point[j]);
		}
		for(int j=1; j<=N; j++){
			knots.emplace_back(point[point.size()-1]+j);
		}
	}
	
	double solve(int n, int i, double x){
    		if (n==0){
    			if(x<=knots[i] || x>knots[i+1]){return 0.0;}
    			else{return 1.0;}
    		}
    		double left = (x-knots[i])/(knots[i+n]-knots[i])*solve(n-1, i, x);
    		double right = (knots[i+n+1]-x)/(knots[i+n+1]-knots[i+1])*solve(n-1, i+1, x);
    		return left+right;
    	}
	
	vector<double> knots;
	int N;
};

class Interpolation_BSpline{
public:
	Interpolation_BSpline(vector<double> point_x, vector<double> point_y, int n) : n(n){
		for(int i=n-1; i>=1; i--){
			knots_x.emplace_back(point_x[0]-i);
		}
		for(int i=0; i<point_x.size(); i++){
			knots_x.emplace_back(point_x[i]);
			knots_y.emplace_back(point_y[i]);
		}
	}
	
	Interpolation_BSpline(vector<double> point, int n, vector<double> coef) : n(n) , coef(coef){
		for(int i=n-1; i>=1; i--){
			knots_x.emplace_back(point[0]-i);
		}
		for(int i=0; i<point.size(); i++){
			knots_x.emplace_back(point[i]);
		}
	}
	
	void natural_BSpline_3_2(){
		int N = knots_x.size();
		B_spline b(knots_x,3);
		MatrixXf A = MatrixXf::Zero(N-2,N-2);
		VectorXf B = VectorXf::Zero(N-2);
		VectorXf temp = VectorXf::Zero(N-2);
		for(int i=0; i<N-2; i++){
			B(i) = knots_y[i];
		}
		for(int i=0; i<N-2; i++){
			for(int j=0; j<N-2; j++){
				if(i==0){
					A(i,i) = (1+(knots_x[3]-knots_x[0])/(knots_x[4]-knots_x[1]))*b.solve(3,0,knots_x[2])+b.solve(3,1,knots_x[2]);
					A(i,i+1) = b.solve(3,2,knots_x[2])-(knots_x[3]-knots_x[0])*(b.solve(3,0,knots_x[2]))/(knots_x[4]-knots_x[1]);
				}else if(i==N-3){
					A(i,i-1) = b.solve(3,N-3,knots_x[N-1])-(b.knots[N+2]-b.knots[N-1])/(b.knots[N+1]-b.knots[N-2]);
					A(i,i) = b.solve(3,N-3,knots_x[N-1])*(1+(b.knots[N+2]-b.knots[N-1])/(b.knots[N+1]-b.knots[N-2]))+b.solve(3,N-2,knots_x[N-1]);
				}else{
					A(i,i-1) = b.solve(3,i,knots_x[i+2]);
					A(i,i) = b.solve(3,i+1,knots_x[i+2]);
					A(i,i+1) = b.solve(3,i+2,knots_x[i+2]);
				}
			}
		}
		temp = A.colPivHouseholderQr().solve(B);
		coef.emplace_back(-(knots_x[3]-knots_x[0])*(temp[1]-temp[0])/(knots_x[4]-knots_x[1])+temp[0]);
		for(int i=0; i<N-2; i++){
			coef.emplace_back(temp[i]);
		}
		coef.emplace_back(temp[N-3]+(b.knots[N+2]-b.knots[N-1])*(temp[N-3]-temp[N-4])/(b.knots[N+1]-b.knots[N-2]));
	}
	
	void complete_Spline_3_2(double boundary_left, double boundary_right){
		int N = knots_x.size();
		B_spline b(knots_x,3);
		MatrixXf A = MatrixXf::Zero(N-2,N-2);
		VectorXf B = VectorXf::Zero(N-2);
		VectorXf temp = VectorXf::Zero(N-2);
		for(int i=0; i<N-2; i++){
			if(i==0){
				B(i) = knots_y[0]+boundary_left*(knots_x[3]-knots_x[0])*b.solve(3,0,knots_x[2])/3/b.solve(2,1,knots_x[2]);
			}else if(i==N-3){
				B(i) = knots_y[N-3]-boundary_right*(b.knots[N+2]-    b.knots[N-1])*b.solve(3,N-1,knots_x[N-1])/3/b.solve(2,N-1,knots_x[N-1]);
			}else{
				B(i) = knots_y[i];
			}
		}
		for(int i=0; i<N-2; i++){
			for(int j=0; j<N-2; j++){
				if(i==0){
					A(i,i) = (b.solve(3,0,knots_x[2])-b.solve(3,0,knots_x[2])*b.solve(2,2,knots_x[2])*(knots_x[3]-knots_x[0])/b.solve(2,1,knots_x[2])/(knots_x[4]-knots_x[1])+b.solve(3,1,knots_x[2]));
					A(i,i+1) = (b.solve(2,2,knots_x[2])*(knots_x[3]-knots_x[0])*b.solve(3,0,knots_x[2])/(knots_x[4]-knots_x[1])/b.solve(2,1,knots_x[2])+b.solve(3,2,knots_x[2]));
				}else if(i==N-3){
					A(i,i-1) = b.solve(3,N-3,knots_x[N-1])+b.solve(2,N-2,knots_x[N-1])*(b.knots[N+2]-b.knots[N-1])*b.solve(3,N-1,knots_x[N-1])/(b.knots[N+1]-b.knots[N-2])/b.solve(2,N-1,knots_x[N-1]);
					A(i,i) = b.solve(3,N-1,knots_x[N-1])+b.solve(3,N-2,knots_x[N-1])-b.solve(2,N-2,knots_x[N-1])*(b.knots[N+2]-b.knots[N-1])*b.solve(3,N-1,knots_x[N-1])/(b.knots[N+1]-b.knots[N-2])/b.solve(2,N-1,knots_x[N-1]);
				}else{
					A(i,i-1) = b.solve(3,i,knots_x[i+2]);
					A(i,i) = b.solve(3,i+1,knots_x[i+2]);
					A(i,i+1) = b.solve(3,i+2,knots_x[i+2]);
				}
			}
		}
		temp = A.colPivHouseholderQr().solve(B);
		coef.emplace_back(temp[0]-(boundary_left+3*(temp[0]-temp[1])*b.solve(2,2,knots_x[2])/(knots_x[4]-knots_x[1]))*(knots_x[3]-knots_x[0])/(3*(b.solve(2,1,knots_x[2]))));
		for(int i=0; i<N-2; i++){
			coef.emplace_back(temp[i]);
		}
		coef.emplace_back(temp[N-3]+(boundary_right-3*(temp[N-3]-temp[N-4])*b.solve(2,N-2,knots_x[N-1])/(b.knots[N+1]-b.knots[N-2]))*(b.knots[N+2]-b.knots[N-1])/(3*(b.solve(2,N-1,knots_x[N-1]))));
	}
	
	void periodic_Spline_3_2(){
		int N = knots_x.size();
		B_spline bb(knots_x,3);
		MatrixXf A = MatrixXf::Zero(N-2,N-2);
		VectorXf B = VectorXf::Zero(N-2);
		VectorXf temp = VectorXf::Zero(N-2);
		for(int i=0; i<N-2; i++){
			B(i) = knots_y[i];
		}
		for(int i=0; i<N-2; i++){
			for(int j=0; j<N-2; j++){
				if(i==0){
					double a = 1/(1+bb.solve(2,N-1,knots_x[N-1])*(bb.knots[N+1]-bb.knots[N-1])/(bb.solve(2,1,knots_x[2]))/(knots_x[3]-knots_x[1]));
					double b = (knots_x[3]-knots_x[0])/(knots_x[4]-knots_x[1]);
					double c = a*b*(knots_x[4]-knots_x[1])*(bb.solve(2,N-1,knots_x[N-1])+bb.solve(2,N-2,knots_x[N-1]))/bb.solve(2,1,knots_x[2])/(bb.knots[N+1]-bb.knots[N-2]);
					A(0,0) = (1+b-a*b-a*b*bb.solve(2,2,knots_x[2])/(bb.solve(2,1,knots_x[2])))*bb.solve(3,0,knots_x[2])+bb.solve(3,1,knots_x[2]);
					A(0,1) = (a*b-b+a*b*bb.solve(2,2,knots_x[2])/bb.solve(2,1,knots_x[2]))*bb.solve(3,0,knots_x[2])+bb.solve(3,2,knots_x[2]);
					A(0,N-4) = c*bb.solve(3,0,knots_x[2]);
					A(0,N-3) = -c*bb.solve(3,0,knots_x[2]);
				}else if(i==N-3){
					double a = 1-bb.solve(2,1,knots_x[2])*(knots_x[3]-knots_x[1])/bb.solve(2,N-1,knots_x[N-1])/(bb.knots[N+1]-bb.knots[N-1]);
					double b = bb.solve(2,1,knots_x[2])*(bb.knots[N+2]-bb.knots[N-1])/bb.solve(2,N-1,knots_x[N-1]);
					double c = (bb.knots[N+2]-bb.knots[N-1])*bb.solve(2,2,knots_x[2])/bb.solve(2,N-1,knots_x[N-1])/(knots_x[4]-knots_x[1]);
					double d = bb.solve(2,N-2,knots_x[N-1])*(bb.knots[N+2]-bb.knots[N-1])/bb.solve(2,N-1,knots_x[N-1])/(bb.knots[N+1]-bb.knots[N-2]);
					double e = b*(knots_x[3]-knots_x[1])/a/(bb.knots[N+1]-bb.knots[N-1])/(bb.knots[N+1]-bb.knots[N-2])-d;
					A(N-3,N-4) = -e;
					A(N-3,N-3) = e+1;
					A(N-3,0) = -b/a/(knots_x[4]-knots_x[1])-c;
					A(N-3,1) = b/a/(knots_x[4]-knots_x[1])+c;
				}else{
					A(i,i-1) = bb.solve(3,i,knots_x[i+2]);
					A(i,i) = bb.solve(3,i+1,knots_x[i+2]);
					A(i,i+1) = bb.solve(3,i+2,knots_x[i+2]);
				}
			}
		}
		temp = A.colPivHouseholderQr().solve(B);
		coef.emplace_back(temp[0]*(1+(knots_x[3]-knots_x[0])/(knots_x[4]-knots_x[1])-((knots_x[3]-knots_x[0])/(knots_x[4]-knots_x[1])*1/(1+bb.solve(2,N-1,knots_x[N-1])*(bb.knots[N+1]-bb.knots[N-1])/(bb.solve(2,1,knots_x[2]))/(knots_x[3]-knots_x[1])))*(1+bb.solve(2,2,knots_x[2])/bb.solve(2,1,knots_x[2])))+temp[1]*(((knots_x[3]-knots_x[0])/(knots_x[4]-knots_x[1])*1/(1+bb.solve(2,N-1,knots_x[N-1])*(bb.knots[N+1]-bb.knots[N-1])/(bb.solve(2,1,knots_x[2]))/(knots_x[3]-knots_x[1])))*(1+bb.solve(2,2,knots_x[2])/bb.solve(2,1,knots_x[2]))-(knots_x[3]-knots_x[0])/(knots_x[4]-knots_x[1]))-((1/(1+bb.solve(2,N-1,knots_x[N-1])*(bb.knots[N+1]-bb.knots[N-1])/(bb.solve(2,1,knots_x[2]))/(knots_x[3]-knots_x[1])))*((knots_x[3]-knots_x[0])/(knots_x[4]-knots_x[1]))*(knots_x[4]-knots_x[1])*(bb.solve(2,N-1,knots_x[N-1])+bb.solve(2,N-2,knots_x[N-1]))/bb.solve(2,1,knots_x[2])/(bb.knots[N+1]-bb.knots[N-2]))*(temp[N-3]-temp[N-4]));
		for(int i=0; i<N-2; i++){
			coef.emplace_back(temp[i]);
		}
		double a = 1-bb.solve(2,1,knots_x[2])*(knots_x[3]-knots_x[1])/bb.solve(2,N-1,knots_x[N-1])/(bb.knots[N+1]-bb.knots[N-1]);
		double b = bb.solve(2,1,knots_x[2])*(bb.knots[N+2]-bb.knots[N-1])/bb.solve(2,N-1,knots_x[N-1]);
		double c = (bb.knots[N+2]-bb.knots[N-1])*bb.solve(2,2,knots_x[2])/bb.solve(2,N-1,knots_x[N-1])/(knots_x[4]-knots_x[1]);
		double d = bb.solve(2,N-2,knots_x[N-1])*(bb.knots[N+2]-bb.knots[N-1])/bb.solve(2,N-1,knots_x[N-1])/(bb.knots[N+1]-bb.knots[N-2]);
		double e = b*(knots_x[3]-knots_x[1])/a/(bb.knots[N+1]-bb.knots[N-1])/(bb.knots[N+1]-bb.knots[N-2])-d;
		coef.emplace_back((e+1)*temp[N-3]-e*temp[N-4]+(c+b/a/(knots_x[4]-knots_x[1]))*temp[1]+(-c-b/a/(knots_x[4]-knots_x[1])))*temp[0];
	}
	
	void BSpline_1_0(){
		for(int i=0; i<knots_x.size(); i++){
			coef.emplace_back(knots_y[i]);
		}
	}
	
	double solve(double x){
		double sum = 0;
		B_spline b(knots_x,n);
		for(int i=0; i<knots_x.size(); i++){
			sum += coef[i]*b.solve(n,i,x);
		}
		return sum;
	}
	
	vector<double> coef;
	
protected:
	int n;
	vector<double> knots_x;
	vector<double> knots_y;
};

double PP_1_0(const vector<double>& knots_x, const vector<double>& knots_y, double& x){
	if(x<knots_x[0] || x>knots_x[knots_x.size()-1]){
		cout << "PP_1_0 :error" << endl;
		return -1;
	}
	int i=0;
	for(i; i<knots_x.size(); i++){
		if(knots_x[i]<=x && x<=knots_x[i+1]){break;}
	}
	return ((knots_y[i+1]-knots_y[i])/(knots_x[i+1]-knots_x[i]))*(x-knots_x[i])+knots_y[i];
};

class PP_3_2_complete{
protected:
	vector<double> x,y;
	vector<double> a,b,c,d;
	VectorXf m;
	double boundary_left;
	double boundary_right;
public:
	PP_3_2_complete(const vector<double>& x, const vector<double>& y, double left, double right) : x(x) , y(y) ,boundary_left(left) , boundary_right(right){
		compute_coefficients();
	}
	
	void compute_coefficients(){
		int n = x.size();
		m.resize(n);
		a.resize(n-1);
		b.resize(n-1);
		c.resize(n-1);
		d.resize(n-1);
		MatrixXf A = MatrixXf::Zero(n-2,n-2);
		VectorXf B = VectorXf::Zero(n-2); 
		VectorXf temp = VectorXf::Zero(n-2);
		for(int i=0; i<n-2; i++){
			if(i!=0){A(i,i-1)=(x[i+2]-x[i+1])/(x[i+2]-x[i]);}
			if(i!=n-3){A(i,i+1)=(x[i+1]-x[i])/(x[i+2]-x[i]);}
			A(i,i) = 2;
		}
		for(int i=0; i<n-2; i++){
			if(i==0){
				B(i) = 3*(x[i+1]-x[i])/(x[i+2]-x[i])*(y[i+2]-y[i+1])/(x[i+2]-x[i])+(x[i+2]-x[i+1])/(x[i+2]-x[i])*(3*(y[i+1]-y[i])/(x[i+1]-x[i])-boundary_left);
				continue;
			}
			if(i==n-3){
				B(i) = (x[i+1]-x[i])/(x[i+2]-x[i])*(3*(y[i+2]-y[i+1])/(x[i+2]-x[i])-boundary_right)+3*(x[i+2]-x[i+1])/(x[i+2]-x[i])*(y[i+1]-y[i])/(x[i+1]-x[i]); 
				continue;
			}
			B(i) = 3*(x[i+1]-x[i])/(x[i+2]-x[i])*(y[i+2]-y[i+1])/(x[i+2]-x[i])+3*(x[i+2]-x[i+1])/(x[i+2]-x[i])*(y[i+1]-y[i])/(x[i+1]-x[i]);
		}
		temp = A.colPivHouseholderQr().solve(B);
		for(int i=0; i<n; i++){
			if(i==0){m[i]=boundary_left;continue;}
			if(i==n-1){m[i]=boundary_right;continue;}
			m[i] = temp[i-1];
		}
		for(int i=0; i<n-1; i++){
			a[i] = y[i];
			b[i] = m[i];
			c[i] = (3*(y[i+1]-y[i])/(x[i+1]-x[i])-2*m[i]-m[i+1])/(x[i+1]-x[i]);
			d[i] = (m[i]+m[i+1]-2*(y[i+1]-y[i])/(x[i+1]-x[i]))/((x[i+1]-x[i])*(x[i+1]-x[i]));
		}
	}
	
	double solve(const double& t){
		if(t<x[0] || t>x[x.size()-1]){cout << "PP_3_2_complete :error" << endl;}
		int i = 0;
		for(i; i<x.size()-1; i++){
			if(x[i]<=t && t<=x[i+1]){break;}
		}
		return a[i]+b[i]*(t-x[i])+c[i]*pow(t-x[i],2)+d[i]*pow(t-x[i],3);
	}
};

class PP_3_2_periodic{
protected:
	vector<double> x,y;
	vector<double> a,b,c,d;
	VectorXf m;
public:
	PP_3_2_periodic(const vector<double>& x, const vector<double>& y) : x(x) , y(y) {
		compute_coefficients(); // y[y.size()-1] == y[0]
	}
	
	void compute_coefficients(){
		int n = x.size();
		m.resize(n);
		a.resize(n-1);
		b.resize(n-1);
		c.resize(n-1);
		d.resize(n-1);
		MatrixXf A = MatrixXf::Zero(n-2,n-2);
		VectorXf B = VectorXf::Zero(n-2);
		VectorXf temp = VectorXf::Zero(n-2);
		for(int i=0; i<n-2; i++){
			if(i!=0){A(i,i-1)=(x[i+2]-x[i+1])/(x[i+2]-x[i]);}
			if(i!=n-3){A(i,i+1)=(x[i+1]-x[i])/(x[i+2]-x[i]);}
			if(i==0){
				A(i,i)=2-(x[i+2]-x[i+1])*(x[n-1]-x[n-2])/(x[i+2]-x[i])/(4*(x[i+1]-x[i])+2*(x[n-1]-x[n-2])); 
				A(i,n-3) = -2*(x[i+2]-x[i+1])*(x[i+1]-x[i])/((x[i+2]-x[i])*(4*(x[i+1]-x[i])+2*(x[n-1]-x[n-2])));
				continue;
			}
			if(i==n-3){
				A(i,i)=2-(x[i+1]-x[i])*(x[1]-x[0])/(x[i+2]-x[i])/(4*(x[1]-x[0])+2*(x[n-1]-x[n-2]));
				A(i,0) = -(x[i+1]-x[i])*(x[n-1]-x[n])/(x[i+2]-x[i])/(4*(x[1]-x[0])+2*(x[n-1]-x[n-2]));
				continue;
			}
			A(i,i) = 2;
		}
		for(int i=0; i<n-2; i++){
			if(i==0){
				B(i) = 3*(x[i+1]-x[i])/(x[i+2]-x[i])*(y[i+2]-y[i+1])/(x[i+2]-x[i])+3*(x[i+2]-x[i+1])/(x[i+2]-x[i])*(y[i+1]-y[i])/(x[i+1]-x[i])-3*(x[i+2]-x[i+1])*(x[n-1]-x[n-2])/(x[i+2]-x[i])*(y[1]-y[0])/(x[1]-x[0])/(4*(x[i+1]-x[i])+2*(x[n-1]-x[n-2]))-6*(x[i+2]-x[i])*(x[i+1]-x[i])/(x[i+2]-x[i])*(y[n-1]-y[n-2])/(x[n-1]-x[n-2])/(4*(x[i+1]-x[i])+2*(x[n-1]-x[n-2]));
				continue;
			}
			if(i==n-3){
				B(i) = 3*(x[i+1]-x[i])/(x[i+2]-x[i])*(y[i+2]-y[i+1])/(x[i+2]-x[i])+3*(x[i+2]-x[i+1])/(x[i+2]-x[i])*(y[i+1]-y[i])/(x[i+1]-x[i])-3*(x[n-2]-x[n-3])*(x[n-1]-x[n-2])*(y[1]-y[0])/(x[n-1]-x[n-3])/(x[1]-x[0])/(4*(x[1]-x[0])+2*(x[n-1]-x[n-2]))-6*(x[n-2]-x[n-3])*(x[1]-x[0])*(y[n-1]-y[n-2])/(x[n-1]-x[n-3])/(x[n-1]-x[n-2])/(4*(x[1]-x[0])+2*(x[n-1]-x[n-2])); 
				continue;
			}
			B(i) = 3*(x[i+1]-x[i])/(x[i+2]-x[i])*(y[i+2]-y[i+1])/(x[i+2]-x[i])+3*(x[i+2]-x[i+1])/(x[i+2]-x[i])*(y[i+1]-y[i])/(x[i+1]-x[i]);
		}
		temp = A.colPivHouseholderQr().solve(B);
		for(int i=0; i<n; i++){
			if(i==0 || i==n-1){
				m[i]=(3*(x[n-1]-x[n-2])*(y[1]-y[0])/(x[1]-x[0])-(x[n-1]-x[n-2])*temp[0]+6*(x[1]-x[0])*(y[n-1]-y[n-2])/(x[n-1]-x[n-2])-2*(x[1]-x[0])*temp[n-3])/(4*(x[1]-x[0])+2*(x[n-1]-x[n-2]));
				continue;
			}
			m[i] = temp[i];
		}
		for(int i=0; i<n-1; i++){
			a[i] = y[i];
			b[i] = m[i];
			c[i] = (3*(y[i+1]-y[i])/(x[i+1]-x[i])-2*m[i]-m[i+1])/(x[i+1]-x[i]);
			d[i] = (m[i]+m[i+1]-2*(y[i+1]-y[i])/(x[i+1]-x[i]))/((x[i+1]-x[i])*(x[i+1]-x[i]));
		}
	}
	
	double solve(const double& t){
		if(t<x[0] || t>x[x.size()-1]){cout << "PP_3_2_periodic :error" << endl;}
		int i = 0;
		for(i; i<x.size()-1; i++){
			if(x[i]<=t && t<=x[i+1]){break;}
		}
		return a[i]+b[i]*(t-x[i])+c[i]*pow(t-x[i],2)+d[i]*pow(t-x[i],3);
	}
};

class PP_3_2_natural{
protected:
	vector<double> x,y;
	vector<double> a,b,c,d;
	VectorXf m;
public:
	PP_3_2_natural(const vector<double>& x, const vector<double>& y) : x(x) , y(y) {
		compute_coefficients();
	}
	
	void compute_coefficients(){
		int n = x.size();
		m.resize(n);
		a.resize(n-1);
		b.resize(n-1);
		c.resize(n-1);
		d.resize(n-1);
		MatrixXf A = MatrixXf::Zero(n-2,n-2);
		VectorXf B = VectorXf::Zero(n-2); 
		VectorXf temp = VectorXf::Zero(n-2);
		for(int i=0; i<n-2; i++){
			if(i!=0){A(i,i-1)=(x[i+2]-x[i+1])/(x[i+2]-x[i]);}
			if(i!=n-3){A(i,i+1)=(x[i+1]-x[i])/(x[i+2]-x[i]);}
			if(i==0){
				A(i,i)=2-(x[i+2]-x[i+1])/(2*(x[i+2]-x[i])); 
				continue;
			}
			if(i==n-3){
				A(i,i)=2-(x[i+1]-x[i])/(2*(x[i+2]-x[i]));
				continue;
			}
			A(i,i) = 2;
		}
		for(int i=0; i<n-2; i++){
			if(i==0){
				B(i) = 3*(x[i+1]-x[i])/(x[i+2]-x[i])*(y[i+2]-y[i+1])/(x[i+2]-x[i])+3/2*(x[i+2]-x[i+1])/(x[i+2]-x[i])*(y[i+1]-y[i])/(x[i+1]-x[i]);
				continue;
			}
			if(i==n-3){
				B(i) = 3/2*(x[i+1]-x[i])/(x[i+2]-x[i])*(y[i+2]-y[i+1])/(x[i+2]-x[i])+3*(x[i+2]-x[i+1])/(x[i+2]-x[i])*(y[i+1]-y[i])/(x[i+1]-x[i]); 
				continue;
			}
			B(i) = 3*(x[i+1]-x[i])/(x[i+2]-x[i])*(y[i+2]-y[i+1])/(x[i+2]-x[i])+3*(x[i+2]-x[i+1])/(x[i+2]-x[i])*(y[i+1]-y[i])/(x[i+1]-x[i]);
		}
		temp = A.colPivHouseholderQr().solve(B);
		for(int i=0; i<n; i++){
			if(i==0){m[i]=(3*(y[1]-y[0])/(x[1]-x[0])-temp[0])/2;continue;}
			if(i==n-1){m[i]=(3*(y[n-1]-y[n-2])/(x[n-1]-x[n-2])-temp[n-3])/2;continue;}
			m[i] = temp[i-1];
		}
		for(int i=0; i<n-1; i++){
			a[i] = y[i];
			b[i] = m[i];
			c[i] = (3*(y[i+1]-y[i])/(x[i+1]-x[i])-2*m[i]-m[i+1])/(x[i+1]-x[i]);
			d[i] = (m[i]+m[i+1]-2*(y[i+1]-y[i])/(x[i+1]-x[i]))/((x[i+1]-x[i])*(x[i+1]-x[i]));
		}
	}
	
	double solve(const double& t){
		if(t<x[0] || t>x[x.size()-1]){cout << "PP_3_2_natural :error" << endl;}
		int i = 0;
		for(i; i<x.size()-1; i++){
			if(x[i]<=t && t<=x[i+1]){break;}
		}
		return a[i]+b[i]*(t-x[i])+c[i]*pow(t-x[i],2)+d[i]*pow(t-x[i],3);
	}
};


