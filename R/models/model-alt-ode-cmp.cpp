#include <R.h>
#include <iostream>

static double parms[11];

#define N parms[0]
#define a parms[1]
#define Tinf parms[2]
#define ldts parms[3]
#define ldte parms[4]
#define betaIn0 parms[5]
#define betaIs0 parms[6]
#define Tin0 parms[7]
#define betaInt parms[8]
#define betaIst parms[9]
#define Tint parms[10]

#define S y[0]
#define E y[1]
#define In y[2]
#define Is y[3]
#define R y[4]

/* initializer */
extern "C" {
  void initmod(void (* odeparms)(int *, double *))
  {
    int c=11;
    odeparms(&c, parms);
  }
}

static double interpolate(double t, double v0, double vt)
{
  if (t < ldts)
    return v0;
  else if (t > ldte)
    return vt;
  else
    return v0 + (t - ldts) / (ldte - ldts) * (vt - v0);
}

/* Derivatives and 1 output variable */
extern "C" {
  void derivs (int *neq, double *t, double *y, double *ydot,
	       double *yout, int *ip)
  {
    if (ip[0] < 2) error("nout should be at least 2");

    if (S < 0 || E < 0 || In < 0 || Is < 0 || R < 0 ||
	S >= N || E >= N || In >= N || Is >= N || R >= N) {
      //std::cerr << S << " " << E << " " << In << " " << Is << " " << R
      //          << std::endl;
      /* dS  */ ydot[0] = 0;
      /* dE  */ ydot[1] = 0;
      /* dIn */ ydot[2] = 0;
      /* dIs */ ydot[3] = 0;
      /* dR  */ ydot[4] = 0;
      return;
    }
   
    const double betaIn = interpolate(*t, betaIn0, betaInt);
    const double betaIs = interpolate(*t, betaIs0, betaIst);
    const double Tin = interpolate(*t, Tin0, Tint);
    const double Tis = Tinf - Tin;
    const double gamman = 1/Tin;
    const double gammas = 1/Tis;
  
    const double inf1 = (betaIn * In) / N * S;
    const double inf2 = (betaIs * Is) / N * S;
    const double got_infected = inf1 + inf2;
    const double got_infectious = a * E;
    const double got_isolated = gamman * In;
    const double got_removed = gammas * Is;
  
    /* dS  */ ydot[0] = -got_infected;
    /* dE  */ ydot[1] = got_infected - got_infectious;
    /* dIn */ ydot[2] = got_infectious - got_isolated;
    /* dIs */ ydot[3] = got_isolated - got_removed;
    /* dR  */ ydot[4] = got_removed;

    yout[0] = Tin * inf1 / In + Tis * inf2 / Is;
    yout[1] = yout[0] * N / S;

  }

  /* The Jacobian matrix */
  void jac(int *neq, double *t, double *y, int *ml, int *mu,
	   double *pd, int *nrowpd, double *yout, int *ip)
  {
    if (S < 0 || E < 0 || In < 0 || Is < 0 || R < 0 ||
	S >= N || E >= N || In >= N || Is >= N || R >= N) {
      //std::cerr << S << " " << E << " " << In << " " << Is << " " << R
      //	  << std::endl;
      return;
    }
   
    const double betaIn = interpolate(*t, betaIn0, betaInt);
    const double betaIs = interpolate(*t, betaIs0, betaIst);
    const double Tin = interpolate(*t, Tin0, Tint);
    const double Tis = Tinf - Tin;
    const double gamman = 1/Tin;
    const double gammas = 1/Tis;

    /*  dS/dS */ pd[               0] = - (betaIn * In + betaIs * Is) / N;
    /*  dE/dS */ pd[               1] = - pd[0];
    /* dIn/dS */ pd[               2] = 0;
    /* dIs/dS */ pd[               3] = 0;
    /*  dR/dS */ pd[               4] = 0;

    /*  dS/dE */  pd[  (*nrowpd)    ] = 0;
    /*  dE/dE */  pd[  (*nrowpd) + 1] = -a;
    /* dIn/dE */  pd[  (*nrowpd) + 2] = a;
    /* dIs/dE */  pd[  (*nrowpd) + 3] = 0;
    /*  dR/dE */  pd[  (*nrowpd) + 4] = 0;

    /*  dS/dIn */ pd[2*(*nrowpd)    ] = betaIn * S / N;
    /*  dE/dIn */ pd[2*(*nrowpd) + 1] = 0;
    /* dIn/dIn */ pd[2*(*nrowpd) + 2] = -gamman;
    /* dIs/dIn */ pd[2*(*nrowpd) + 3] = gamman;
    /*  dR/dIn */ pd[2*(*nrowpd) + 4] = 0;

    /*  dS/dIs */ pd[3*(*nrowpd)    ] = betaIs * S / N;
    /*  dE/dIs */ pd[3*(*nrowpd) + 1] = 0;
    /* dIn/dIs */ pd[3*(*nrowpd) + 2] = 0;
    /* dIs/dIs */ pd[3*(*nrowpd) + 3] = -gammas;
    /*  dR/dIs */ pd[3*(*nrowpd) + 4] = gammas;

    /*  dS/dR */  pd[4*(*nrowpd)    ] = 0;
    /*  dE/dR */  pd[4*(*nrowpd) + 1] = 0;
    /* dIn/dR */  pd[4*(*nrowpd) + 2] = 0;
    /* dIs/dR */  pd[4*(*nrowpd) + 3] = 0;
    /*  dR/dR */  pd[4*(*nrowpd) + 4] = 0;
  }
}
