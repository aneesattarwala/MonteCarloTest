
//  main.cpp
//  D.a
//  Copyright (c) 2015 MichaelScott. All rights reserved.
// We wish to add functionality to the Monte Carlo pricer by providing estimates for the standard deviation (SD) and standard error (SE)
// TESTMC.cpp

#include "OptionData.hpp"
#include "UtilitiesDJD/RNG/NormalGenerator.hpp"
#include "UtilitiesDJD/Geometry/Range.cpp"
#include "UtilitiesDJD/VectorsAndMatrices/ArrayMechanisms.cpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
using namespace std;
// Function to calculate standard deviation and standard error of MC simulation
tuple<double, double> AdvMonteCarlo (const vector<double>& P, double T, double r)
{
    double SD, SE, _Sum=0, _SqrSum=0;  // Standard deviation and standard error
    double M = P.size();
    vector<double>::const_iterator it;
    for (it=P.begin();it!=P.end();++it)
    {
        _Sum +=*it;
        _SqrSum += (*it)*(*it);
    }
    SD = sqrt((_SqrSum - _Sum*_Sum/M))/(M-1)*exp(-2*r*T);
    SE = SD/sqrt(M);
    return  tuple <double, double> (SD,SE);
    
}

template <class T> void print(const std::vector<T>& myList)
{  // A generic print function for vectors
	
	std::cout << std::endl << "Size of vector is " << myList.size() << "\n[";
    
	// We must use a const iterator here, otherwise we get a compiler error.
	typename std::vector<T>::const_iterator i;
	for (i = myList.begin(); i != myList.end(); ++i)
	{
        std::cout << *i << ",";
        
	}
    
	std::cout << "]\n";
}

namespace SDEDefinition
{ // Defines drift + diffusion + data
    
	OptionData* data;				// The data for the option MC
    
	double drift(double t, double X)
	{ // Drift term
        
		return (data->r)*X; // r - D
	}
    
	
	double diffusion(double t, double X)
	{ // Diffusion term
        
		double betaCEV = 1.0;
		return data->sig * pow(X, betaCEV);
		
	}
    
	double diffusionDerivative(double t, double X)
	{ // Diffusion term, needed for the Milstein method
        
		double betaCEV = 1.0;
		return 0.5 * (data->sig) * (betaCEV) * pow(X, 2.0 * betaCEV - 1.0);
	}
} // End of namespace


int main()
{
	std::cout <<  "1 factor MC with explicit Euler\n";
	OptionData myOption;                   // Batch 1 parameters
//	myOption.K = 65.0;
//	myOption.T = 0.25;
//	myOption.r = 0.08;
//	myOption.sig = 0.3;
//	myOption.type = -1;	// Put -1, Call +1
//	double S_0 = 60;
//	
    
    myOption.K = 100.0;                // batch 2 parameters
	myOption.T = 1;
	myOption.r = 0.0;
	myOption.sig = 0.2;
	myOption.type = 1;	// Put -1, Call +1
	double S_0 = 100;
    
    
	long N = 100;
	std::cout << "Number of subintervals in time: ";
	std::cin >> N;
    
	// Create the basic SDE (Context class)
	Range<double> range (0.0, myOption.T);
	double VOld = S_0;
	double VNew;
    
	Vector<double, long> x = range.mesh(N);
	
    
	// V2 mediator stuff
	long NSim = 50000;
	std::cout << "Number of simulations: ";
	std::cin >> NSim;
    
	double k = myOption.T / double (N);
	double sqrk = sqrt(k);
    
	// Normal random number
	double dW;
	double price = 0.0;	// Option price
    vector<double> payoff;
	// NormalGenerator is a base class
	NormalGenerator* myNormal = new BoostNormal();
    
	using namespace SDEDefinition;
	data = &myOption;
    
	Vector<double> res;
	int coun = 0; // Number of times S hits origin
    
	// A.
	for (long i = 1; i <= NSim; ++i)
	{ // Calculate a path at each iteration
        
//		if ((i/10000) * 10000 == i)
//		{// Give status after each 1000th iteration
//            
//            std::cout << i << std::endl;
//		}
        
		VOld = S_0;
		for (long index = x.MinIndex()+1; index <= x.MaxIndex(); ++index)
		{
            
			// Create a random number
			dW = myNormal->getNormal();
            
			// The FDM (in this case explicit Euler)
			VNew = VOld  + (k * drift(x[index-1], VOld))
            + (sqrk * diffusion(x[index-1], VOld) * dW);
            
			VOld = VNew;
            
			// Spurious values
			if (VNew <= 0.0) coun++;
		}
        
		double tmp = myOption.myPayOffFunction(VNew);
		payoff.push_back(tmp);
        price += (tmp)/double(NSim);
	}
	
	// D. Finally, discounting the average price
	price *= exp(-myOption.r * myOption.T);
    
	// Cleanup; V2 use scoped pointer
	//delete myNormal;
    
	cout << "Price, after discounting: " << price << ", " << std::endl;
	cout << "Number of times origin is hit: " << coun << endl;
    
    tuple<double,double> St_Result;    // tuple for storing simulation statistics
    
    St_Result = AdvMonteCarlo(payoff, myOption.T, myOption.r);  // function call
    
    cout << "Standard Deviation :"<< get<0>(St_Result)<<endl;
    cout << "Standard Error :"<< get<1>(St_Result)<<endl;

	return 0;
}