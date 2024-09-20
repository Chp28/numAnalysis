#include "Function.hpp"
#include "EquationSolver.hpp"
#include <iostream>
#include <cmath>
using namespace std;

class F : public Function{
public:
	double operator() (double x) const{
		return tan(x);
	}
};

void solve(){
	cout << "Solving tan(x) on [4.5,7.7]" <<endl;
	Newton_Method solver(F(),7);
	double x = solver.solve();
	cout << "A root is: " << x <<endl;
	return;
}

int main(){
	solve();
}
