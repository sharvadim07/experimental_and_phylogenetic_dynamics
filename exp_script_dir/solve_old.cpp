#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>


#include <iostream>
#include <fstream>
#include <utility>

#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

// global num of equations
//#define dim 12
// global num of parameters
//#define param_num 15
//#define alpha 1

double a0[param_num + 1];
double a[param_num + 1];


#define Min(a, b) ((a < b) ? a : b)
#define Max(a, b) ((a > b) ? a : b)


struct f {
	void operator()(const vector_type &x , vector_type &res , double t) {
		res[0] = 0;
		// derivs Equations start 
		//
		// derivs Equations end 
		
		// res[1] = a[1]*x[1]*(1 - (x[1] + x[2] + x[3] + x[4])/a[2]) - a[3]*x[1]*x[5] - a[4]*x[1]*x[6] - a[5]*x[1];
		// res[2] = a[3]*x[1]*x[5] - a[4]*x[2]*x[6] - a[6]*x[2];
		// res[3] = a[4]*x[1]*x[6] - a[3]*x[3]*x[5] - a[6]*x[3];
		// res[4] = a[4]*x[2]*x[6] + a[3]*x[3]*x[5] - a[6]*x[4];
		// res[5] = - a[7]*x[5]*(x[1] + x[2] + x[3] + x[4]) + a[9]*Max((x[7]*x[11])/(x[2] + alpha), 0) + a[9]*Max((x[9]*x[12])/(x[4] + alpha) , 0); 
		// res[6] = - a[8]*x[6]*(x[1] + x[2] + x[3] + x[4]) + a[10]*Max((x[10]*x[12])/(x[4] + alpha), 0);  
		// res[7] = a[3]*x[5]*(x[1] + x[2]) + a[11]*x[7]*Max(Min((1 - x[7]/(a[13]*x[2] + alpha)), 1), 0) - a[9]*Max((x[7]*x[11])/(x[2] + alpha), 0) - a[4]*x[7]*x[6] - a[6]*x[7]; 
		// res[8] = a[4]*x[6]*(x[1] + x[3]) + a[12]*x[8]*Max(Min((1 - x[8]/(a[13]*x[3] + alpha)), 1), 0) - a[3]*x[8]*x[5] - a[6]*x[8]; 
		// res[9] = a[3]*x[5]*(x[3] + x[4]) + a[11]*x[9]*Max(Min((1 - (x[9] + x[10])/(a[13]*x[4] + alpha)), 1), 0) - a[9]*Max((x[9]*x[12])/(x[4] + alpha), 0) + a[4]*x[7]*x[6] - a[6]*x[9]; 
		// res[10] = a[4]*x[6]*(x[2] + x[4]) + a[12]*x[10]*Max(Min((1 - (x[9] + x[10])/(a[13]*x[4] + alpha)), 1), 0) - a[10]*Max((x[10]*x[12])/(x[4] + alpha), 0) + a[3]*x[8]*x[5] - a[6]*x[10]; 
		// res[11] = a[14]*x[7]*Max(Min((1 - x[11]/(a[15]*x[2] + alpha)), 1), 0) - a[9]*Max((x[7]*x[11])/(x[2] + alpha), 0) - a[4]*x[11]*x[6] - a[6]*x[11]; 
		// res[12] = a[14]*x[9]*Max(Min((1 - x[12]/(a[15]*x[4] + alpha)), 1), 0) - a[9]*Max((x[9]*x[12])/(x[4] + alpha), 0) - a[10]*Max((x[10]*x[12])/(x[4] + alpha), 0) + a[4]*x[11]*x[6] - a[6]*x[12];  
	}
};

struct f_jacobi {
	void operator()(const vector_type &x, matrix_type &J, const double &t, vector_type &dfdt) {
		
		//jac Jacobian start
		//jac Jacobian end
		
		// J(0, 0) = 0;
		// J(0, 1) = 0;
		// J(0, 2) = 0;
		// J(0, 3) = 0;
		// J(0, 4) = 0;
		// J(0, 5) = 0;
		// J(0, 6) = 0;
		// J(0, 7) = 0;
		// J(0, 8) = 0;
		// J(0, 9) = 0;
		// J(0, 10) = 0;
		// J(0, 11) = 0;
		// J(0, 12) = 0;

		// J(1, 0) = 0;
		// J(1, 1) = a[1]*(1 - (x[1] + x[2] + x[3] + x[4])/a[2]) - a[1]*x[1]/a[2] - a[3]*x[5] - a[4]*x[6] - a[5];
		// J(1, 2) = -a[1]*x[1]/a[2];
		// J(1, 3) = -a[1]*x[1]/a[2];
		// J(1, 4) = -a[1]*x[1]/a[2];
		// J(1, 5) = -a[3]*x[1];
		// J(1, 6) = -a[4]*x[1];
		// J(1, 7) = 0;
		// J(1, 8) = 0;
		// J(1, 9) = 0;
		// J(1, 10) = 0;
		// J(1, 11) = 0;
		// J(1, 12) = 0;

		// J(2, 0) = 0;
		// J(2, 1) = a[3]*x[5];
		// J(2, 2) = -a[4]*x[6] - a[6];
		// J(2, 3) = 0;
		// J(2, 4) = 0;
		// J(2, 5) = a[3]*x[1];
		// J(2, 6) = -a[4]*x[2];
		// J(2, 7) = 0;
		// J(2, 8) = 0;
		// J(2, 9) = 0;
		// J(2, 10) = 0;
		// J(2, 11) = 0;
		// J(2, 12) = 0;

		// J(3, 0) = 0;
		// J(3, 1) = a[4]*x[6];
		// J(3, 2) = 0;
		// J(3, 3) = -a[3]*x[5] - a[6];
		// J(3, 4) = 0;
		// J(3, 5) = -a[3]*x[3];
		// J(3, 6) = a[4]*x[1];
		// J(3, 7) = 0;
		// J(3, 8) = 0;
		// J(3, 9) = 0;
		// J(3, 10) = 0;
		// J(3, 11) = 0;
		// J(3, 12) = 0;

		// J(4, 0) = 0;
		// J(4, 1) = 0;
		// J(4, 2) = a[4]*x[6];
		// J(4, 3) = a[3]*x[5];
		// J(4, 4) = -a[6];
		// J(4, 5) = a[3]*x[3];
		// J(4, 6) = a[4]*x[2];
		// J(4, 7) = 0;
		// J(4, 8) = 0;
		// J(4, 9) = 0;
		// J(4, 10) = 0;
		// J(4, 11) = 0;
		// J(4, 12) = 0;

		// J(5, 0) = 0;
		// J(5, 1) = a[7]*x[5];
		// J(5, 2) = a[7]*x[5] - a[9]*x[11]*x[7]/pow(alpha + x[2], 2);
		// J(5, 3) = a[7]*x[5];
		// J(5, 4) = a[7]*x[5] - a[9]*x[12]*x[9]/pow(alpha + x[4], 2);
		// J(5, 5) = a[7]*(x[1] + x[2] + x[3] + x[4]);
		// J(5, 6) = 0;
		// J(5, 7) = a[9]*x[11]/(alpha + x[2]);
		// J(5, 8) = 0;
		// J(5, 9) = a[9]*x[12]/(alpha + x[4]);
		// J(5, 10) = 0;
		// J(5, 11) = a[9]*x[7]/(alpha + x[2]);
		// J(5, 12) = a[9]*x[9]/(alpha + x[4]);

		// J(6, 0) = 0;
		// J(6, 1) = a[8]*x[6];
		// J(6, 2) = a[8]*x[6];
		// J(6, 3) = a[8]*x[6];
		// J(6, 4) = -a[10]*x[10]*x[12]/pow(alpha + x[4], 2) + a[8]*x[6];
		// J(6, 5) = 0;
		// J(6, 6) = a[8]*(x[1] + x[2] + x[3] + x[4]);
		// J(6, 7) = 0;
		// J(6, 8) = 0;
		// J(6, 9) = 0;
		// J(6, 10) = a[10]*x[12]/(alpha + x[4]);
		// J(6, 11) = 0;
		// J(6, 12) = a[10]*x[10]/(alpha + x[4]);

		// J(7, 0) = 0;
		// J(7, 1) = a[3]*x[5];
		// J(7, 2) = a[11]*a[13]*pow(x[7], 2)/pow(a[13]*x[2] + alpha, 2) + a[3]*x[5] + a[9]*x[11]*x[7]/pow(alpha + x[2], 2);
		// J(7, 3) = 0;
		// J(7, 4) = 0;
		// J(7, 5) = a[3]*(x[1] + x[2]);
		// J(7, 6) = -a[4]*x[7];
		// J(7, 7) = -a[11]*x[7]/(a[13]*x[2] + alpha) + a[11]*(-x[7]/(a[13]*x[2] + alpha) + 1) - a[4]*x[6] - a[6] - a[9]*x[11]/(alpha + x[2]);
		// J(7, 8) = 0;
		// J(7, 9) = 0;
		// J(7, 10) = 0;
		// J(7, 11) = -a[9]*x[7]/(alpha + x[2]);
		// J(7, 12) = 0;

		// J(8, 0) = 0;
		// J(8, 1) = a[4]*x[6];
		// J(8, 2) = 0;
		// J(8, 3) = a[12]*a[13]*pow(x[8], 2)/pow(a[13]*x[3] + alpha, 2) + a[4]*x[6];
		// J(8, 4) = 0;
		// J(8, 5) = -a[3]*x[8];
		// J(8, 6) = a[4]*(x[1] + x[3]);
		// J(8, 7) = 0;
		// J(8, 8) = -a[12]*x[8]/(a[13]*x[3] + alpha) + a[12]*(-x[8]/(a[13]*x[3] + alpha) + 1) - a[3]*x[5] - a[6];
		// J(8, 9) = 0;
		// J(8, 10) = 0;
		// J(8, 11) = 0;
		// J(8, 12) = 0;

		// J(9, 0) = 0;
		// J(9, 1) = 0;
		// J(9, 2) = 0;
		// J(9, 3) = a[3]*x[5];
		// J(9, 4) = a[11]*a[13]*x[9]*(x[10] + x[9])/pow(a[13]*x[4] + alpha, 2) + a[3]*x[5] + a[9]*x[12]*x[9]/pow(alpha + x[4], 2);
		// J(9, 5) = a[3]*(x[3] + x[4]);
		// J(9, 6) = a[4]*x[7];
		// J(9, 7) = a[4]*x[6];
		// J(9, 8) = 0;
		// J(9, 9) = -a[11]*x[9]/(a[13]*x[4] + alpha) + a[11]*(-(x[10] + x[9])/(a[13]*x[4] + alpha) + 1) - a[6] - a[9]*x[12]/(alpha + x[4]);
		// J(9, 10) = -a[11]*x[9]/(a[13]*x[4] + alpha);
		// J(9, 11) = 0;
		// J(9, 12) = -a[9]*x[9]/(alpha + x[4]);

		// J(10, 0) = 0;
		// J(10, 1) = 0;
		// J(10, 2) = a[4]*x[6];
		// J(10, 3) = 0;
		// J(10, 4) = a[10]*x[10]*x[12]/pow(alpha + x[4], 2) + a[12]*a[13]*x[10]*(x[10] + x[9])/pow(a[13]*x[4] + alpha, 2) + a[4]*x[6];
		// J(10, 5) = a[3]*x[8];
		// J(10, 6) = a[4]*(x[2] + x[4]);
		// J(10, 7) = 0;
		// J(10, 8) = a[3]*x[5];
		// J(10, 9) = -a[12]*x[10]/(a[13]*x[4] + alpha);
		// J(10, 10) = -a[10]*x[12]/(alpha + x[4]) - a[12]*x[10]/(a[13]*x[4] + alpha) + a[12]*(-(x[10] + x[9])/(a[13]*x[4] + alpha) + 1) - a[6];
		// J(10, 11) = 0;
		// J(10, 12) = -a[10]*x[10]/(alpha + x[4]);

		// J(11, 0) = 0;
		// J(11, 1) = 0;
		// J(11, 2) = a[14]*a[15]*x[11]*x[7]/pow(a[15]*x[2] + alpha, 2) + a[9]*x[11]*x[7]/pow(alpha + x[2], 2);
		// J(11, 3) = 0;
		// J(11, 4) = 0;
		// J(11, 5) = 0;
		// J(11, 6) = -a[4]*x[11];
		// J(11, 7) = a[14]*(-x[11]/(a[15]*x[2] + alpha) + 1) - a[9]*x[11]/(alpha + x[2]);
		// J(11, 8) = 0;
		// J(11, 9) = 0;
		// J(11, 10) = 0;
		// J(11, 11) = -a[14]*x[7]/(a[15]*x[2] + alpha) - a[4]*x[6] - a[6] - a[9]*x[7]/(alpha + x[2]);
		// J(11, 12) = 0;

		// J(12, 0) = 0;
		// J(12, 1) = 0;
		// J(12, 2) = 0;
		// J(12, 3) = 0;
		// J(12, 4) = a[10]*x[10]*x[12]/pow(alpha + x[4], 2) + a[14]*a[15]*x[12]*x[9]/pow(a[15]*x[4] + alpha, 2) + a[9]*x[12]*x[9]/pow(alpha + x[4], 2);
		// J(12, 5) = 0;
		// J(12, 6) = a[4]*x[11];
		// J(12, 7) = 0;
		// J(12, 8) = 0;
		// J(12, 9) = a[14]*(-x[12]/(a[15]*x[4] + alpha) + 1) - a[9]*x[12]/(alpha + x[4]);
		// J(12, 10) = -a[10]*x[12]/(alpha + x[4]);
		// J(12, 11) = a[4]*x[6];
		// J(12, 12) = -a[10]*x[10]/(alpha + x[4]) - a[14]*x[9]/(a[15]*x[4] + alpha) - a[6] - a[9]*x[9]/(alpha + x[4]);
		
		//dfdt start
		//dfdt end
		// dfdt[0] = 0;
		// dfdt[1] = 0;
		// dfdt[2] = 0;
		// dfdt[3] = 0;
		// dfdt[4] = 0;
		// dfdt[5] = 0;
		// dfdt[6] = 0;
		// dfdt[7] = 0;
		// dfdt[8] = 0;
		// dfdt[9] = 0;
		// dfdt[10] = 0;
		// dfdt[11] = 0;
		// dfdt[12] = 0;
		
	}
};

void print_iter(const vector_type &x , const double t) {
	printf("%f\t", t);
	for (int i = 1; i < dim + 1; i++)
		cout << abs(x[i]) << '\t';
	
	if (x[1] + x[2] + x[3] + x[4] < 1)
		exit(0);
	
	cout << endl;
}


pid_t child_pid = -1;

void catchAlarm(int sig) {
    //std::cerr << "You will now die, Mr Solver!!\n";
    kill(child_pid,SIGKILL);
}

int main(int argc, char** argv) {

	// Die in 30 seconds....
        //signal(SIGALRM, catchAlarm);
        //alarm(30); 

	struct sigaction sact;
	sigemptyset(&sact.sa_mask);
	sact.sa_flags = 0;
	sact.sa_handler = catchAlarm;
	sigaction(SIGALRM, &sact, NULL);
	//alarm(30);
	

	if (argc != param_num + 5) {
		cerr << argc << endl;
		cerr << "Please provide parameters (t1, t2, step, accuracy, a1, a2, ...)" << endl;
		return -1;
	}

	// Solving ODE on interval [t1, t2] with step h
	double t1 = atof(argv[1]);
	double t2 = atof(argv[2]);
	double h = atof(argv[3]);
	double accuracy = atof(argv[4]);
	
	a0[0] = 0;
	a[0] = 0;
	for (int i = 1; i < param_num + 1; i++) {
		a0[i] = atof(argv[4 + i]);
		a[i] = a0[i];
	}

	// print header traj file
	//cout << "#time\tHC\tWC\tDC\tWDC\tWV\tDV\tLW_W\tLD_D\tLW_WD\tLD_WD\tP_W\tP_WD" << endl;
	
	vector_type x(dim + 1);	
	
	x[0] = 0.; 
	//initial variables values start
	//initial variables values end
	// x[1] = 40000.;
	// x[2] = 0.;
	// x[3] = 0.; 
	// x[4] = 0.; 
	// x[5] = 500.;
	// x[6] = 5000.;
	// x[7] = 0.;
	// x[8] = 0.;
	// x[9] = 0.;
	// x[10] = 0.;
	// x[11] = 0.;
	// x[12] = 0.;
	    //clock_t tStart = clock();
	    child_pid = fork();
	    int st = 0;
	    if (child_pid == 0) {
	            integrate_const(make_dense_output< rosenbrock4< double > >( accuracy , accuracy ), 
	              make_pair(f(), f_jacobi()), 
	              x, t1, t2, h, print_iter);
	              exit(1);
		}
	    else
		{
		alarm(60);
		wait(&st);
		}
	return 0;
}
