#include <iostream>
#include <R.h>

static const int PARAM_N = 60;
static double parms[PARAM_N];

#define Ny         parms[0]
#define No         parms[1]
#define a1         parms[2]
#define a2         parms[3]
#define gamma      parms[4]
#define eta        parms[5]
#define tuncertain parms[6]
#define funcertain parms[7]
#define t0         parms[8] // Start of physical distancing beta0
#define t1         parms[9] // Start of lockdown            beta1
#define t2         parms[10] // End of lockdown transition   beta2
#define betay0     parms[11]
#define betao0     parms[12]
#define betayo0    parms[13]
#define betay1     parms[14]
#define betao1     parms[15]
#define betayo1    parms[16]
#define betay2     parms[17]
#define betao2     parms[18]
#define betayo2    parms[19]
#define t3         parms[20]      // linear change point
#define betay3     parms[21]
#define betao3     parms[22]
#define betayo3    parms[23]
#define t4         parms[24]      // linear change point
#define betay4     parms[25]
#define betao4     parms[26]
#define betayo4    parms[27]
#define t5         parms[28]      // step change point, lockdown 2
#define betay5     parms[29]
#define betao5     parms[30]
#define betayo5    parms[31]
#define t6         parms[32]      // step change point, end of lockdown 2
#define betay6     parms[33]
#define betao6     parms[34]
#define betayo6    parms[35]
#define t7         parms[36]      // linear change point (sept 15?)
#define betay7     parms[37]
#define betao7     parms[38]
#define betayo7    parms[39]
#define t8         parms[40]      // linear change point, oct9 rules
#define betay8     parms[41]
#define betao8     parms[42]
#define betayo8    parms[43]
#define t9         parms[44]      // linear change point, oct16 rules
#define betay9     parms[45]
#define betao9     parms[46]
#define betayo9    parms[47]
#define t10        parms[48]      // linear change point
#define betay10    parms[49]
#define betao10    parms[50]
#define betayo10   parms[51]
#define t11        parms[52]      // linear change point
#define betay11    parms[53]
#define betao11    parms[54]
#define betayo11   parms[55]
#define t12        parms[56]      // step change point
#define betay12    parms[57]
#define betao12    parms[58]
#define betayo12   parms[59]

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
    // std::cerr << "=== initmod ===" << std::endl;
  }
}

static double interpolate(double t,
			  double vt0, double vt1, double vt2, double vt3,
			  double vt4, double vt5, double vt6, double vt7,
			  double vt8, double vt9, double vt10, double vt11,
			  double vt12)
{
  if (t > t12)
    return vt12;
  else if (t > t11)
    return vt11;
  else if (t > t10)
    return vt10 + (t - t10) / (t11 - t10) * (vt11 - vt10);
  else if (t > t9)
    return vt9;
  else if (t > t8)
    return vt8 + (t - t8) / (t9 - t8) * (vt9 - vt8);
  else if (t > t7)
    return vt7;
  else if (t > t6)
    return vt6;
  else if (t > t5)
    return vt5;
  else if (t > t4)
    return vt4;
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

static double uncertain(double t, double v)
{
  if (t > tuncertain)
    return v * funcertain;
  else
    return v;
}

/* Derivatives and 1 output variable */
extern "C" {
  void derivs (int *neq, double *t, double *y, double *ydot,
	       double *yout, int *ip)
  {
    if (ip[0] < 9) error("nout should be at least 9");
    
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

    // std::cerr << "derivs t=" << *t << std::endl;

    const double mbetay = interpolate(*t,
				      betay0, betay1, betay2, betay3, betay4,
				      betay5, betay6, betay7, betay8, betay9,
				      betay10, betay11, betay12);
    const double mbetao = interpolate(*t,
				      betao0, betao1, betao2, betao3, betao4,
				      betao5, betao6, betao7, betao8, betao9,
				      betao10, betao11, betao12);
    const double mbetayo = interpolate(*t,
				       betayo0, betayo1, betayo2, betayo3, betayo4,
				       betayo5, betayo6, betayo7, betayo8, betayo9,
				       betayo10, betayo11, betayo12);

    const double betay = uncertain(*t, mbetay);
    const double betao = uncertain(*t, mbetao);
    const double betayo = uncertain(*t, mbetayo);
   
    const double Syeff = std::max(0.0, Ny * 0.8 - (Ny - Sy));

    const double ygot_infected = (betay * Iy) / Ny * Syeff;
    const double ygot_latent2 = a1 * E1y;
    const double ygot_infectious = a2 * E2y;
    const double ygot_removed = gamma * Iy;
    const double ygot_reverted = eta * Ry;

    const double ogot_infected = ((betao * Io) / No + (betayo * Iy) / Ny) * So;
    const double ogot_latent2 = a1 * E1o;
    const double ogot_infectious = a2 * E2o;
    const double ogot_removed = gamma * Io;
    const double ogot_reverted = eta * Ro;

    /* dSy  */ ydot[0] = -ygot_infected + ygot_reverted;
    /* dE1y */ ydot[1] = ygot_infected - ygot_latent2;
    /* dE2y */ ydot[2] = ygot_latent2 - ygot_infectious;
    /* dIy  */ ydot[3] = ygot_infectious - ygot_removed;
    /* dRy  */ ydot[4] = ygot_removed - ygot_reverted;

    /* dSo  */ ydot[5] = -ogot_infected + ogot_reverted;
    /* dE1o */ ydot[6] = ogot_infected - ogot_latent2;
    /* dE2o */ ydot[7] = ogot_latent2 - ogot_infectious;
    /* dIo  */ ydot[8] = ogot_infectious - ogot_removed;
    /* dRo  */ ydot[9] = ogot_removed - ogot_reverted;

    yout[0] = 1/gamma * (ygot_infected + ogot_infected) / (Iy + Io);
    yout[1] = yout[0] * (Ny + No) / (Sy + So);

    yout[2] = 1/gamma * ((betay * Iy) / Ny * Sy + (betayo * Iy) / Ny * So) / Iy;
    yout[3] = 1/gamma *                            (betao * Io) / No * So / Io;

    yout[4] = ygot_infected;
    yout[5] = ogot_infected;
    yout[6] = betay;
    yout[7] = betao;
    yout[8] = betayo;
  }
}
