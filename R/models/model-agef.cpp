#include <iostream>
#include <R.h>

static const int PARAM_N = 37;
static double parms[PARAM_N];

#define Ny parms[0]
#define No parms[1]
#define a1 parms[2]
#define a2 parms[3]
#define gamma parms[4]
#define t0 parms[5] // Start of physical distancing beta0
#define t1 parms[6] // Start of lockdown            beta1
#define t2 parms[7] // End of lockdown transition   beta2
#define betay0 parms[8]
#define betao0 parms[9]
#define betayo0 parms[10]
#define betay1 parms[11]
#define betao1 parms[12]
#define betayo1 parms[13]
#define betay2 parms[14]
#define betao2 parms[15]
#define betayo2 parms[16]
#define t3 parms[17]      // linear change point
#define betay3 parms[18]
#define betao3 parms[19]
#define betayo3 parms[20]
#define t4 parms[21]      // linear change point
#define betay4 parms[22]
#define betao4 parms[23]
#define betayo4 parms[24]
#define t5 parms[25]      // linear change point, lockdown 2
#define betay5 parms[26]
#define betao5 parms[27]
#define betayo5 parms[28]
#define t6 parms[29]      // linear change point, end of lockdown 2 transition
#define betay6 parms[30]
#define betao6 parms[31]
#define betayo6 parms[32]
#define t7 parms[33]      // step transition
#define betay7 parms[34]
#define betao7 parms[35]
#define betayo7 parms[36]

#define Sy y[0]
#define E1y y[1]
#define E2y y[2]
#define Iy y[3]
#define Ry y[4]
#define So y[5]
#define E1o y[6]
#define E2o y[7]
#define Io y[8]
#define Ro y[9]

/* initializer */
extern "C" {
  void initmod(void (* odeparms)(int *, double *))
  {
    int c=PARAM_N;
    odeparms(&c, parms);
  }
}

static double interpolate(double t,
			  double vt0, double vt1, double vt2, double vt3,
			  double vt4, double vt5, double vt6, double vt7)
{
  if (t > t7)
    return vt7;
  if (t > t6)
    return vt6;
  else if (t > t5)
    return vt5 + (t - t5) / (t6 - t5) * (vt6 - vt5);
  else if (t > t4)
    return vt4 + (t - t4) / (t5 - t4) * (vt5 - vt4);
  else if (t > t3)
    return vt3 + (t - t3) / (t4 - t3) * (vt4 - vt3);
  else if (t > t2)
    return vt2 + (t - t2) / (t3 - t2) * (vt3 - vt2);
  else if (t > t1)
    return vt1 + (t - t1) / (t2 - t1) * (vt2 - vt1);
  else if (t > t0)
    return vt0 + (t - t0) / (t1 - t0) * (vt1 - vt0);
  else
    return vt0;
}

/* Derivatives and 1 output variable */
extern "C" {
  void derivs (int *neq, double *t, double *y, double *ydot,
	       double *yout, int *ip)
  {
    if (ip[0] < 4) error("nout should be at least 4");
    
    if (Sy < 0  || E1y < 0  || E2y < 0  || Iy < 0  || Ry < 0  ||
	Sy > Ny || E1y > Ny || E2y > Ny || Iy > Ny || Ry > Ny ||
	So < 0  || E1o < 0  || E2o < 0  || Io < 0  || Ro < 0  ||
	So > No || E1o > No || E2o > No || Io > No || Ro > No) {
      /* dSy  */ ydot[0] = 0;
      /* dE1y */ ydot[1] = 0;
      /* dE2y */ ydot[2] = 0;
      /* dIy  */ ydot[3] = 0;
      /* dRy  */ ydot[4] = 0;
      /* dSo  */ ydot[5] = 0;
      /* dE1o */ ydot[6] = 0;
      /* dE2o */ ydot[7] = 0;
      /* dIo  */ ydot[8] = 0;
      /* dRo  */ ydot[9] = 0;
      return;
    }

    const double betay = interpolate(*t,
				     betay0, betay1, betay2, betay3, betay4,
				     betay5, betay6, betay7);
    const double betao = interpolate(*t,
				     betao0, betao1, betao2, betao3, betao4,
				     betao5, betao6, betao7);
    const double betayo = interpolate(*t,
				      betayo0, betayo1, betayo2, betayo3, betayo4,
				      betayo5, betayo6, betayo7);

    const double ygot_infected = (betay * Iy) / Ny * Sy;
    const double ygot_latent2 = a1 * E1y;
    const double ygot_infectious = a2 * E2y;
    const double ygot_removed = gamma * Iy;

    const double ogot_infected = ((betao * Io) / No + (betayo * Iy) / Ny) * So;
    const double ogot_latent2 = a1 * E1o;
    const double ogot_infectious = a2 * E2o;
    const double ogot_removed = gamma * Io;

    /* dSy  */ ydot[0] = -ygot_infected;
    /* dE1y */ ydot[1] = ygot_infected - ygot_latent2;
    /* dE2y */ ydot[2] = ygot_latent2 - ygot_infectious;
    /* dIy  */ ydot[3] = ygot_infectious - ygot_removed;
    /* dRy  */ ydot[4] = ygot_removed;

    /* dSo  */ ydot[5] = -ogot_infected;
    /* dE1o */ ydot[6] = ogot_infected - ogot_latent2;
    /* dE2o */ ydot[7] = ogot_latent2 - ogot_infectious;
    /* dIo  */ ydot[8] = ogot_infectious - ogot_removed;
    /* dRo  */ ydot[9] = ogot_removed;

    yout[0] = 1/gamma * (ygot_infected + ogot_infected) / (Iy + Io);
    yout[1] = yout[0] * (Ny + No) / (Sy + So);

    yout[2] = 1/gamma * ((betay * Iy) / Ny * Sy + (betayo * Iy) / Ny * So) / Iy;
    yout[3] = 1/gamma *                            (betao * Io) / No * So / Io;
  }
}
