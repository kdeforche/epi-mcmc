#include <R.h>

static double parms[8];

#define N parms[0]
#define a parms[1]
#define ldts parms[2]
#define ldte parms[3]
#define beta0 parms[4]
#define Tinf0 parms[5]
#define betat parms[6]
#define Tinft parms[7]

#define S y[0]
#define E y[1]
#define I y[2]
#define R y[3]

/* initializer */
extern "C" {
  void initmod(void (* odeparms)(int *, double *))
  {
    int c=8;
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
    
    if (S < 0 || E < 0 || I < 0 || R < 0 ||
	S >= N || E >= N || I >= N || R >= N) {
      //std::cerr << S << " " << E << " " << In << " " << Is << " " << R
      //          << std::endl;
      /* dS  */ ydot[0] = 0;
      /* dE  */ ydot[1] = 0;
      /* dI  */ ydot[2] = 0;
      /* dI  */ ydot[3] = 0;
      return;
    }

    const double beta = interpolate(*t, beta0, betat);
    const double Tinf = interpolate(*t, Tinf0, Tinft);
    const double gamma = 1/Tinf;
  
    const double inf = (beta * I) / N * S;
    const double got_infected = inf;
    const double got_infectious = a * E;
    const double got_removed = gamma * I;
  
    /* dS  */ ydot[0] = -got_infected;
    /* dE  */ ydot[1] = got_infected - got_infectious;
    /* dI  */ ydot[2] = got_infectious - got_removed;
    /* dR  */ ydot[3] = got_removed;

    yout[0] = Tinf * inf / I;
    yout[1] = yout[0] * N / S;
  }

  /* The Jacobian matrix */
  void jac(int *neq, double *t, double *y, int *ml, int *mu,
	   double *pd, int *nrowpd, double *yout, int *ip)
  {
    if (S < 0 || E < 0 || I < 0 || R < 0 ||
	S >= N || E >= N || I >= N || R >= N) {
      //std::cerr << S << " " << E << " " << In << " " << Is << " " << R
      //	  << std::endl;
      return;
    }

    const double beta = interpolate(*t, beta0, betat);
    const double gamma = 1/interpolate(*t, Tinf0, Tinft);
  
    /*  dS/dS */ pd[               0] = -beta * I / N;
    /*  dE/dS */ pd[               1] = -pd[0];
    /*  dI/dS */ pd[               2] = 0;
    /*  dR/dS */ pd[               3] = 0;

    /*  dS/dE */ pd[  (*nrowpd)    ] = 0;
    /*  dE/dE */ pd[  (*nrowpd) + 1] = -a;
    /*  dI/dE */ pd[  (*nrowpd) + 2] = a;
    /*  dR/dE */ pd[  (*nrowpd) + 3] = 0;

    /*  dS/dI */ pd[2*(*nrowpd)    ] = beta * S / N;
    /*  dE/dI */ pd[2*(*nrowpd) + 1] = 0;
    /*  dI/dI */ pd[2*(*nrowpd) + 2] = -gamma;
    /*  dR/dI */ pd[2*(*nrowpd) + 3] = gamma;

    /*  dS/dR */ pd[3*(*nrowpd)    ] = 0;
    /*  dE/dR */ pd[3*(*nrowpd) + 1] = 0;
    /*  dI/dR */ pd[3*(*nrowpd) + 2] = 0;
    /*  dR/dR */ pd[3*(*nrowpd) + 3] = 0;
  }
}
