#include <R.h>

static double parms[11];

#define N parms[0]
#define a parms[1]
#define phts parms[2]
#define ldts parms[3]
#define ldte parms[4]
#define betat0 parms[5]
#define Tinft0 parms[6]
#define betat1 parms[7]
#define Tinft1 parms[8]
#define betat2 parms[9]
#define Tinft2 parms[10]

#define S y[0]
#define E y[1]
#define I y[2]
#define R y[3]

/* initializer */
extern "C" {
  void initmod(void (* odeparms)(int *, double *))
  {
    int c=11;
    odeparms(&c, parms);
  }
}

static double interpolate(double t, double vt0, double vt1, double vt2)
{
  /* If phts > ldts: make vt0 ~ vt1 */
  
  if (t < phts)
    return vt0;
  else if (t > ldts) {
    if (t > ldte)
      return vt2;
    else
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

    const double beta = interpolate(*t, betat0, betat1, betat2);
    const double Tinf = interpolate(*t, Tinft0, Tinft1, Tinft2);
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

    const double beta = interpolate(*t, betat0, betat1, betat2);
    const double gamma = 1/interpolate(*t, Tinft0, Tinft1, Tinft2);
  
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
