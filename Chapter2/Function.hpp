#ifndef FUNCTION
#define FUNCTION

class Function{
public:
	virtual double operator() (double x) const = 0;
	virtual double derivative(double x) const {
		return (operator()(x+1e-6)-operator()(x))/1e-6;
	}
};

#endif
