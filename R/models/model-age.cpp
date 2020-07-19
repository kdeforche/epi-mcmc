#include <R.h>

static const int PARAM_N = 24;
static double parms[PARAM_N];

#define Ny parms[0]
#define No parms[1]
#define a parms[2]
#define gamma parms[3]
#define phts parms[4] // Start of physical distancing beta0
#define ldts parms[5] // Start of lockdown            beta1
#define ldte parms[6] // End of lockdown transition   beta2
#define betay0 parms[7]
#define betao0 parms[8]
#define betayo0 parms[9]
#define betay1 parms[10]
#define betao1 parms[11]
#define betayo1 parms[12]
#define betay2 parms[13]
#define betao2 parms[14]
#define betayo2 parms[15]
#define t3 parms[16]
#define betay3 parms[17]
#define betao3 parms[18]
#define betayo3 parms[19]
#define t4 parms[20]
#define betay4 parms[21]
#define betao4 parms[22]
#define betayo4 parms[23]

#define Sy y[0]
#define Ey y[1]
#define Iy y[2]
#define Ry y[3]
#define So y[4]
#define Eo y[5]
#define Io y[6]
#define Ro y[7]

/* initializer */
extern "C" {
  void initmod(void (* odeparms)(int *, double *))
  {
    int c=PARAM_N;
    odeparms(&c, parms);
  }
}

static double interpolate(double t, double vt0, double vt1, double vt2, double vt3, double vt4)
{
  if (t < phts)
    return vt0;
  else if (t > ldts) {
    if (t > ldte) {
      if (t > t3)
	if (t > t4)
	  return vt4;
	else
	  return vt3;
      else
	return vt2;
    } else
      return vt1 + (t - ldts) / (ldte - ldts) * (vt2 - vt1);
  } else {
    return vt0 + (t - phts) / (ldts - phts) * (vt1 - vt0);
  }
}


/* Derivatives and 1 output variable */
extern "C" {
  void derivs (int *neq, double *t, double *y, double *ydot,
	       double *yout, int *ip)
  {
    if (ip[0] < 2) error("nout should be at least 2");
    
    if (Sy < 0   || Ey < 0   || Iy < 0   || Ry < 0   ||
	Sy >= Ny || Ey >= Ny || Iy >= Ny || Ry >= Ny ||
	So < 0   || Eo < 0   || Io < 0   || Ro < 0   ||
	So >= No || Eo >= No || Io >= No || Ro >= No) {
      /* dSy  */ ydot[0] = 0;
      /* dEy  */ ydot[1] = 0;
      /* dIy  */ ydot[2] = 0;
      /* dIy  */ ydot[3] = 0;
      /* dSo  */ ydot[4] = 0;
      /* dEo  */ ydot[5] = 0;
      /* dIo  */ ydot[6] = 0;
      /* dIo  */ ydot[7] = 0;
      return;
    }

    const double betay = interpolate(*t, betay0, betay1, betay2, betay3, betay4);
    const double betao = interpolate(*t, betao0, betao1, betao2, betao3, betao4);
    const double betayo = interpolate(*t, betayo0, betayo1, betayo2, betayo3, betayo4);
  
    const double ygot_infected = (betay * Iy) / Ny * Sy;
    const double ygot_infectious = a * Ey;
    const double ygot_removed = gamma * Iy;

    const double ogot_infected = ((betao * Io) / No + (betayo * Iy) / Ny) * So;
    const double ogot_infectious = a * Eo;
    const double ogot_removed = gamma * Io;

    /* dSy  */ ydot[0] = -ygot_infected;
    /* dEy  */ ydot[1] = ygot_infected - ygot_infectious;
    /* dIy  */ ydot[2] = ygot_infectious - ygot_removed;
    /* dRy  */ ydot[3] = ygot_removed;

    /* dSo  */ ydot[4] = -ogot_infected;
    /* dEo  */ ydot[5] = ogot_infected - ogot_infectious;
    /* dIo  */ ydot[6] = ogot_infectious - ogot_removed;
    /* dRo  */ ydot[7] = ogot_removed;

    yout[0] = 1/gamma * (ygot_infected + ogot_infected) / (Iy + Io);
    yout[1] = yout[0] * (Ny + No) / (Sy + So);
  }
}
