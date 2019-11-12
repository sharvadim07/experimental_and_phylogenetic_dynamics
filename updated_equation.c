#include <R.h>
#include <math.h>
#include <stdio.h>
#define Min(a, b) ((a < b) ? a : b)
#define Max(a, b) ((a > b) ? a : b)
#define J(a, b) pd[a*(*nrowpd) + b]

static double a_log[17];
double alpha = 1;

// global num of parameters
#define NUM_OF_PAR 17

/* initializer */
void initmod(void (* odeparms)(int *, double *)) {    
	int N = NUM_OF_PAR;
    odeparms(&N, a_log);
}
/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *x, double *res, double *yout, int *ip) {
    //if (ip[0] < 1) error("nout should be at least 1");
    double a[NUM_OF_PAR];
    for (int i = 0; i < NUM_OF_PAR; i++)
        a[i] = pow(10, a_log[i]);

    res[0] = 0;
    // derivs Equations start 
res[1]=a[1]*x[1]-a[1]*x[1]*(x[1]+x[2]+x[3])/(x[1]+x[2]+x[3]+a[14])-a[2]*x[1]-a[3]*x[4]*x[1]-a[4]*x[5]*x[1];
res[2]=a[4]*x[5]*x[1]-a[3]*x[4]*x[2]-a[11]*x[2];
res[3]=a[4]*x[5]*x[3]+a[3]*x[4]*x[2]+a[3]*x[4]*x[1]-a[11]*x[3];
res[4]=a[7]*x[7]*x[9]/x[3]-a[3]*a[5]*x[4]*(x[1]+x[2]+x[3]);
res[5]=a[8]*x[8]*x[9]/x[3]-a[4]*a[6]*x[5]*(x[1]+x[2]+x[3]);
res[6]=a[4]*x[5]*(x[1]+x[2])-a[3]*x[4]*x[6]-a[11]*x[6];
res[7]=a[3]*x[4]*(x[1]+x[2]+x[3])+a[9]*x[7]-(a[9]*x[7]*(x[7]+x[8])/x[3])/((x[7]+x[8])/x[3]+a[15])-a[7]*x[7]*x[9]/x[3]-a[11]*x[7];
res[8]=a[4]*x[5]*x[3]+a[10]*x[8]*x[7]/(x[7]+x[8])-((a[10]*x[8]*x[7])*(x[7]+x[8])/x[3])/((x[7]+x[8])*((x[7]+x[8])/x[3]+a[15]))+a[13]*x[7]-(a[13]*x[7]*(x[7]+x[8])/x[3])/((x[7]+x[8])/x[3]+a[15])+a[3]*x[4]*x[6]-a[8]*x[8]*x[9]/x[3]-a[11]*x[8];
res[9]=a[12]*x[7]-(a[12]*x[7]*x[9]/x[3])/(x[9]/x[3]+a[16])-a[7]*x[7]*x[9]/x[3]-a[8]*x[8]*x[9]/x[3]-a[11]*x[9];
    // derivs Equations end 

	// derivs Sum of vars
yout[0]=x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7]+x[8]+x[9];
}
/* The Jacobian matrix */
void jac(int *neq, double *t, double *x, int *ml, int *mu, double *pd, int *nrowpd, double *yout, int *ip) {
	double a[NUM_OF_PAR];
    for (int i = 0; i < NUM_OF_PAR; i++)
        a[i] = pow(10, a_log[i]);
		
		//jac Jacobian start
J(0, 0) = 0;
J(0, 1) = 0;
J(0, 2) = 0;
J(0, 3) = 0;
J(0, 4) = 0;
J(0, 5) = 0;
J(0, 6) = 0;
J(0, 7) = 0;
J(0, 8) = 0;
J(0, 9) = 0;
J(1, 0) = 0;
J(1, 1) = a[1] - a[2] - x[5]*a[4] - x[4]*a[3] - x[1]*a[1]/(x[1] + x[2] + a[14] + x[3]) - a[1]*(x[1] + x[2] + x[3])/(x[1] + x[2] + a[14] + x[3]) + x[1]*a[1]*(x[1] + x[2] + x[3])/pow((x[1] + x[2] + a[14] + x[3]),2);
J(1, 2) = -x[1]*a[1]/(x[1] + x[2] + a[14] + x[3]) + x[1]*a[1]*(x[1] + x[2] + x[3])/pow((x[1] + x[2] + a[14] + x[3]),2);
J(1, 3) = -x[1]*a[1]/(x[1] + x[2] + a[14] + x[3]) + x[1]*a[1]*(x[1] + x[2] + x[3])/pow((x[1] + x[2] + a[14] + x[3]),2);
J(1, 4) = -x[1]*a[3];
J(1, 5) = -x[1]*a[4];
J(1, 6) = 0;
J(1, 7) = 0;
J(1, 8) = 0;
J(1, 9) = 0;
J(2, 0) = 0;
J(2, 1) = x[5]*a[4];
J(2, 2) = -a[11] - x[4]*a[3];
J(2, 3) = 0;
J(2, 4) = -x[2]*a[3];
J(2, 5) = x[1]*a[4];
J(2, 6) = 0;
J(2, 7) = 0;
J(2, 8) = 0;
J(2, 9) = 0;
J(3, 0) = 0;
J(3, 1) = x[4]*a[3];
J(3, 2) = x[4]*a[3];
J(3, 3) = -a[11] + x[5]*a[4];
J(3, 4) = x[1]*a[3] + x[2]*a[3];
J(3, 5) = a[4]*x[3];
J(3, 6) = 0;
J(3, 7) = 0;
J(3, 8) = 0;
J(3, 9) = 0;
J(4, 0) = 0;
J(4, 1) = -x[4]*a[5]*a[3];
J(4, 2) = -x[4]*a[5]*a[3];
J(4, 3) = -x[4]*a[5]*a[3] - x[7]*x[9]*a[7]/pow(x[3],2);
J(4, 4) = -a[5]*a[3]*(x[1] + x[2] + x[3]);
J(4, 5) = 0;
J(4, 6) = 0;
J(4, 7) = x[9]*a[7]/x[3];
J(4, 8) = 0;
J(4, 9) = x[7]*a[7]/x[3];
J(5, 0) = 0;
J(5, 1) = -x[5]*a[6]*a[4];
J(5, 2) = -x[5]*a[6]*a[4];
J(5, 3) = -x[5]*a[6]*a[4] - x[8]*x[9]*a[8]/pow(x[3],2);
J(5, 4) = 0;
J(5, 5) = -a[6]*a[4]*(x[1] + x[2] + x[3]);
J(5, 6) = 0;
J(5, 7) = 0;
J(5, 8) = x[9]*a[8]/x[3];
J(5, 9) = x[8]*a[8]/x[3];
J(6, 0) = 0;
J(6, 1) = x[5]*a[4];
J(6, 2) = x[5]*a[4];
J(6, 3) = 0;
J(6, 4) = -x[6]*a[3];
J(6, 5) = (x[1] + x[2])*a[4];
J(6, 6) = -a[11] - x[4]*a[3];
J(6, 7) = 0;
J(6, 8) = 0;
J(6, 9) = 0;
J(7, 0) = 0;
J(7, 1) = x[4]*a[3];
J(7, 2) = x[4]*a[3];
J(7, 3) = x[4]*a[3] + x[7]*x[9]*a[7]/pow(x[3],2) + x[7]*(x[7] + x[8])*a[9]/(pow(x[3],2)*(a[15] + (x[7] + x[8])/x[3])) - x[7]*pow((x[7] + x[8]),2)*a[9]/pow((x[3],3)*pow((a[15] + (x[7] + x[8])/x[3]),2));
J(7, 4) = a[3]*(x[1] + x[2] + x[3]);
J(7, 5) = 0;
J(7, 6) = 0;
J(7, 7) = -a[11] + a[9] - x[9]*a[7]/x[3] - x[7]*a[9]/(x[3]*(a[15] + (x[7] + x[8])/x[3])) - (x[7] + x[8])*a[9]/(x[3]*(a[15] + (x[7] + x[8])/x[3])) + x[7]*(x[7] + x[8])*a[9]/pow((x[3],2)*pow((a[15] + (x[7] + x[8])/x[3]),2));
J(7, 8) = -x[7]*a[9]/(x[3]*(a[15] + (x[7] + x[8])/x[3])) + x[7]*(x[7] + x[8])*a[9]/pow((x[3],2)*pow((a[15] + (x[7] + x[8])/x[3]),2));
J(7, 9) = -x[7]*a[7]/x[3];
J(8, 0) = 0;
J(8, 1) = 0;
J(8, 2) = 0;
J(8, 3) = x[5]*a[4] + x[8]*x[9]*a[8]/pow(x[3],2) + x[7]*x[8]*a[10]/(pow(x[3],2)*(a[15] + (x[7] + x[8])/x[3])) + x[7]*(x[7] + x[8])*a[13]/(pow(x[3],2)*(a[15] + (x[7] + x[8])/x[3])) - x[7]*pow((x[7] + x[8]),2)*a[13]/pow((x[3],3)*pow((a[15] + (x[7] + x[8])/x[3]),2)) - x[7]*x[8]*(x[7] + x[8])*a[10]/pow((x[3],3)*pow((a[15] + (x[7] + x[8])/x[3]),2));
J(8, 4) = x[6]*a[3];
J(8, 5) = a[4]*x[3];
J(8, 6) = x[4]*a[3];
J(8, 7) = a[13] + x[8]*a[10]/(x[7] + x[8]) - x[7]*x[8]*a[10]/pow((x[7] + x[8]),2) - x[7]*a[13]/(x[3]*(a[15] + (x[7] + x[8])/x[3])) - x[8]*a[10]/(x[3]*(a[15] + (x[7] + x[8])/x[3])) - (x[7] + x[8])*a[13]/(x[3]*(a[15] + (x[7] + x[8])/x[3])) + x[7]*x[8]*a[10]/pow((x[3],2)*pow((a[15] + (x[7] + x[8])/x[3]),2)) + x[7]*(x[7] + x[8])*a[13]/pow((x[3],2)*pow((a[15] + (x[7] + x[8])/x[3]),2));
J(8, 8) = -a[11] + x[7]*a[10]/(x[7] + x[8]) - x[9]*a[8]/x[3] - x[7]*x[8]*a[10]/pow((x[7] + x[8]),2) - x[7]*a[13]/(x[3]*(a[15] + (x[7] + x[8])/x[3])) - x[7]*a[10]/(x[3]*(a[15] + (x[7] + x[8])/x[3])) + x[7]*x[8]*a[10]/pow((x[3],2)*pow((a[15] + (x[7] + x[8])/x[3]),2)) + x[7]*(x[7] + x[8])*a[13]/pow((x[3],2)*pow((a[15] + (x[7] + x[8])/x[3]),2));
J(8, 9) = -x[8]*a[8]/x[3];
J(9, 0) = 0;
J(9, 1) = 0;
J(9, 2) = 0;
J(9, 3) = x[7]*x[9]*a[7]/pow(x[3],2) + x[8]*x[9]*a[8]/pow(x[3],2) + x[7]*x[9]*a[12]/(pow(x[3],2)*(a[16] + x[9]/x[3])) - x[7]*pow(x[9],2)*a[12]/pow((x[3],3)*pow((a[16] + x[9]/x[3]),2));
J(9, 4) = 0;
J(9, 5) = 0;
J(9, 6) = 0;
J(9, 7) = a[12] - x[9]*a[7]/x[3] - x[9]*a[12]/(x[3]*(a[16] + x[9]/x[3]));
J(9, 8) = -x[9]*a[8]/x[3];
J(9, 9) = -a[11] - x[7]*a[7]/x[3] - x[8]*a[8]/x[3] - x[7]*a[12]/(x[3]*(a[16] + x[9]/x[3])) + x[7]*x[9]*a[12]/pow((x[3],2)*pow((a[16] + x[9]/x[3]),2));
		
		/*
        J(0, 0) = 0;
                J(0, 1) = 0;
                J(0, 2) = 0;
                J(0, 3) = 0;
                J(0, 4) = 0;
                J(0, 5) = 0;
                J(0, 6) = 0;
                J(0, 7) = 0;

		J(1, 0) = 0;
		J(1, 1) = a[1]*(1 - (x[1] + x[2] + x[3] + x[4] + x[7])/a[2]) - a[1]*x[1]/a[2] - a[5]   // x[13] added
                          - (1.0-a[16])*a[1]*x[7]/a[2]  // Added
                          - a[3]*x[5] - a[4]*x[6]
                          - (a[3]*a[14]*x[5] + a[4]*a[15]*x[6]) * ((2*x[1]*(x[1] + x[2] + x[3] + x[4] + x[7] + alpha) - x[1]*x[1]) / (pow (x[1] + x[2] + x[3] + x[4] + x[7] + alpha, 2)))  // Added
                          - (a[12]*a[14]*x[5]*x[7] + a[13]*a[15]*x[6]*x[7]) * (((x[1] + x[2] + x[3] + x[4] + x[7] + alpha) - x[1]) / (pow (x[1] + x[2] + x[3] + x[4] + x[7] + alpha, 2))); // Added
		J(1, 2) = -a[1]*x[1]/a[2]
                          - (1.0-a[16])*a[1]*x[7]/a[2]  // Added
                          + (a[3]*a[14]*x[5] + a[4]*a[15]*x[6]) * pow ((x[1] / (x[1] + x[2] + x[3] + x[4] + x[7] + alpha)), 2)   // Added
                          + (a[12]*a[14]*x[5]*x[7] + a[13]*a[15]*x[6]*x[7]) * x[1] / (pow ((x[1] + x[2] + x[3] + x[4] + x[7] + alpha), 2));  // Added
		J(1, 3) = -a[1]*x[1]/a[2]
                          - (1.0-a[16])*a[1]*x[7]/a[2]  // Added
                          + (a[3]*a[14]*x[5] + a[4]*a[15]*x[6]) * pow ((x[1] / (x[1] + x[2] + x[3] + x[4] + x[7] + alpha)), 2)  //Added
                          + (a[12]*a[14]*x[5]*x[7] + a[13]*a[15]*x[6]*x[7]) * x[1] / (pow ((x[1] + x[2] + x[3] + x[4] + x[7] + alpha), 2));  // Added
		J(1, 4) = -a[1]*x[1]/a[2]
                          - (1.0-a[16])*a[1]*x[7]/a[2]  // Added
                          + (a[3]*a[14]*x[5] + a[4]*a[15]*x[6]) * pow ((x[1] / (x[1] + x[2] + x[3] + x[4] + x[7] + alpha)), 2)   // Added
                          + (a[12]*a[14]*x[5]*x[7] + a[13]*a[15]*x[6]*x[7]) * x[1] / (pow ((x[1] + x[2] + x[3] + x[4] + x[7] + alpha), 2));  // Added
		J(1, 5) = -a[3]*x[1]
                          - (a[3]*a[14]*x[1] + a[12]*a[14]*x[7]) * x[1] / (x[1] + x[2] + x[3] + x[4] + x[7] + alpha);  // Added
		J(1, 6) = -a[4]*x[1]
                          - (a[4]*a[15]*x[1] + a[13]*a[15]*x[7]) * x[1] / (x[1] + x[2] + x[3] + x[4] + x[7] + alpha);  // Added;
                J(1, 7) = -a[1]*x[1]/a[2]    // Added
                           + (1.0-a[16])*a[1]*(1 - (x[1] + x[2] + x[3] + x[4] + x[7])/a[2]) - (1.0-a[16])*a[1]*x[7] / a[2]            // Added
                           + (a[3]*a[14]*x[5] + a[4]*a[15]*x[6]) * x[1]*x[1] / (pow ((x[1] + x[2] + x[3] + x[4] + x[7] + alpha), 2))   // Added
                           + (a[12]*a[14]*x[5] + a[13]*a[15]*x[6]) * x[1] * (((x[1] + x[2] + x[3] + x[4] + x[7] + alpha) - x[7])   / (pow ((x[1] + x[2] + x[3] + x[4] + x[7] + alpha), 2)));  // Added

		J(2, 0) = 0;
		J(2, 1) = a[3]*x[5];
		J(2, 2) = -a[4]*x[6] - a[6];
		J(2, 3) = 0;
		J(2, 4) = 0;
		J(2, 5) = a[3]*x[1] + a[12]*x[7];
		J(2, 6) = -a[4]*x[2];
                J(2, 7) = a[12]*x[5];

		J(3, 0) = 0;
		J(3, 1) = a[4]*x[6];
		J(3, 2) = 0;
		J(3, 3) = -a[3]*x[5] - a[6];
		J(3, 4) = 0;
		J(3, 5) = -a[3]*x[3];
		J(3, 6) = a[4]*x[1] + a[13]*x[7];
                J(3, 7) = a[13]*x[6];

		J(4, 0) = 0;
		J(4, 1) = 0;
		J(4, 2) = a[4]*x[6];
		J(4, 3) = a[3]*x[5];
		J(4, 4) = -a[6];
		J(4, 5) = a[3]*x[3];
		J(4, 6) = a[4]*x[2];
		J(4, 7) = 0;

                J(5, 0) = 0;
                J(5, 1) = -a[7]*x[5];
                J(5, 2) = -a[7]*x[5] + a[9];
                J(5, 3) = -a[7]*x[5];
                J(5, 4) = -a[7]*x[5] + a[11];
                J(5, 5) = -a[7]*(x[1] + x[2] + x[3] + x[4] + x[7]);
                J(5, 6) = 0;
                J(5, 7) = -a[7]*x[5];

                J(6, 0) = 0;
                J(6, 1) = -a[8]*x[6];
                J(6, 2) = -a[8]*x[6];
                J(6, 3) = -a[8]*x[6];
                J(6, 4) = -a[8]*x[6] + a[10];
                J(6, 5) = 0;
                J(6, 6) = -a[8]*(x[1] + x[2] + x[3] + x[4] + x[7]);
                J(6, 7) = -a[8]*x[6];

                J (7, 0) = 0;
	        J (7, 1) = -a[16]*a[1]*x[7]/a[2]
                            + a[14] * (((x[1] + x[2] + x[3] + x[4] + x[7] + alpha) - x[1]) / (pow ((x[1] + x[2] + x[3] + x[4] + x[7] + alpha), 2))) * (a[3]*x[1]*x[5] + a[12]*x[7]*x[5])
                            + a[14] * (x[1] / (x[1] + x[2] + x[3] + x[4] + x[7] + alpha)) * a[3]*x[5]
                            + a[15] * (((x[1] + x[2] + x[3] + x[4] + x[7] + alpha) - x[1]) / (pow ((x[1] + x[2] + x[3] + x[4] + x[7] + alpha), 2))) * (a[4]*x[1]*x[6] + a[13]*x[7]*x[6])
                            + a[15] * (x[1] / (x[1] + x[2] + x[3] + x[4] + x[7] + alpha)) * a[4]*x[6];
                J (7, 2) = -a[16]*a[1]*x[7]/a[2]
                            - a[14] * (x[1] / pow ((x[1] + x[2] + x[3] + x[4] + x[7] + alpha), 2)) * (a[3]*x[1]*x[5] + a[12]*x[7]*x[5])
                            - a[15] * (x[1] / pow ((x[1] + x[2] + x[3] + x[4] + x[7] + alpha), 2)) * (a[4]*x[1]*x[6] + a[13]*x[7]*x[6]);
                J (7, 3) = -a[16]*a[1]*x[7]/a[2]
                            - a[14] * (x[1] / pow ((x[1] + x[2] + x[3] + x[4] + x[7] + alpha), 2)) * (a[3]*x[1]*x[5] + a[12]*x[7]*x[5])
                            - a[15] * (x[1] / pow ((x[1] + x[2] + x[3] + x[4] + x[7] + alpha), 2)) * (a[4]*x[1]*x[6] + a[13]*x[7]*x[6]);
                J (7, 4) = -a[16]*a[1]*x[7]/a[2]
                            - a[14] * (x[1] / pow ((x[1] + x[2] + x[3] + x[4] + x[7] + alpha), 2)) * (a[3]*x[1]*x[5] + a[12]*x[7]*x[5])
                            - a[15] * (x[1] / pow ((x[1] + x[2] + x[3] + x[4] + x[7] + alpha), 2)) * (a[4]*x[1]*x[6] + a[13]*x[7]*x[6]);
                J (7, 5) = -a[12]*x[7]
                            + a[14] * (x[1] / (x[1] + x[2] + x[3] + x[4] + x[7] + alpha)) * (a[3]*x[1] + a[12]*x[7]);
                J (7, 6) = -a[13]*x[7]
                            + a[15] * (x[1] / (x[1] + x[2] + x[3] + x[4] + x[7] + alpha)) * (a[4]*x[1] + a[13]*x[7]);
	        J (7, 7) = a[16]*a[1]*(1 - (x[1] + x[2] + x[3] + x[4] + x[7])/a[2]) - a[16]*a[1]*x[7]/a[2] - a[5]
                             - a[12]*x[5] - a[13]*x[6]
                             - a[14] * (x[1] / pow ((x[1] + x[2] + x[3] + x[4] + x[7] + alpha), 2)) * (a[3]*x[1]*x[5] + a[12]*x[7]*x[5])
                             + a[14] * (x[1] / (x[1] + x[2] + x[3] + x[4] + x[7] + alpha)) * a[12]*x[5]
                             - a[15] * (x[1] / pow ((x[1] + x[2] + x[3] + x[4] + x[7] + alpha), 2)) * (a[4]*x[1]*x[6] + a[13]*x[7]*x[6])
                             + a[15] * (x[1] / (x[1] + x[2] + x[3] + x[4] + x[7] + alpha)) * a[13]*x[6];
		*/
}
