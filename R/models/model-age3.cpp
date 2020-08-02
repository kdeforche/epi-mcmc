#include <R.h>

#include <iostream>

static const int PARAM_N = 45;
static double parms[PARAM_N];

#define Ny parms[0]
#define Nm parms[1]
#define No parms[2]
#define a parms[3]
#define gamma parms[4]
#define t0 parms[5] // Start of lockdown            beta0
#define t1 parms[6] // End of lockdown transition   beta1

#define betay0 parms[7]
#define betam0 parms[8]
#define betao0 parms[9]
#define betaym0 parms[10]
#define betayo0 parms[11]
#define betamy0 parms[12]
#define betamo0 parms[13]

#define betay1 parms[14]
#define betam1 parms[15]
#define betao1 parms[16]
#define betaym1 parms[17]
#define betayo1 parms[18]
#define betamy1 parms[19]
#define betamo1 parms[20]

#define t2 parms[21]
#define betay2 parms[22]
#define betam2 parms[23]
#define betao2 parms[24]
#define betaym2 parms[25]
#define betayo2 parms[26]
#define betamy2 parms[27]
#define betamo2 parms[28]

#define t3 parms[29]
#define betay3 parms[30]
#define betam3 parms[31]
#define betao3 parms[32]
#define betaym3 parms[33]
#define betayo3 parms[34]
#define betamy3 parms[35]
#define betamo3 parms[36]

#define t4 parms[37]
#define betay4 parms[38]
#define betam4 parms[39]
#define betao4 parms[40]
#define betaym4 parms[41]
#define betayo4 parms[42]
#define betamy4 parms[43]
#define betamo4 parms[44]

#define Sy y[0]
#define Ey y[1]
#define Iy y[2]
#define Ry y[3]
#define Sm y[4]
#define Em y[5]
#define Im y[6]
#define Rm y[7]
#define So y[8]
#define Eo y[9]
#define Io y[10]
#define Ro y[11]

/* initializer */
extern "C" {
  void initmod(void (* odeparms)(int *, double *))
  {
    int c=PARAM_N;
    odeparms(&c, parms);
  }
}

static double interpolate(double t,
			  double vt0, double vt1, double vt2, double vt3, double vt4)
{
  if (t > t4)
    return vt4;
  else if (t > t3)
    return vt3;
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
    if (ip[0] < 5) error("nout should be at least 5");
    
    if (Sy < 0   || Ey < 0   || Iy < 0   || Ry < 0   ||
	Sy > Ny || Ey > Ny || Iy > Ny || Ry > Ny ||
	Sm < 0   || Em < 0   || Im < 0   || Rm < 0   ||
	Sm > Nm || Em > Nm || Im > Nm || Rm > Nm ||
	So < 0   || Eo < 0   || Io < 0   || Ro < 0   ||
	So > No || Eo > No || Io > No || Ro > No) {
      std::cerr << "Oops ?" << std::endl;
      /* dSy  */ ydot[0] = 0;
      /* dEy  */ ydot[1] = 0;
      /* dIy  */ ydot[2] = 0;
      /* dIy  */ ydot[3] = 0;
      /* dSm  */ ydot[4] = 0;
      /* dEm  */ ydot[5] = 0;
      /* dIm  */ ydot[6] = 0;
      /* dIm  */ ydot[7] = 0;
      /* dSo  */ ydot[8] = 0;
      /* dEo  */ ydot[9] = 0;
      /* dIo  */ ydot[10] = 0;
      /* dIo  */ ydot[11] = 0;
      return;
    }

    const double betay = interpolate(*t,
				     betay0, betay1, betay2, betay3, betay4);
    const double betam = interpolate(*t,
				     betam0, betam1, betam2, betam3, betam4);
    const double betao = interpolate(*t,
				     betao0, betao1, betao2, betao3, betao4);
    const double betaym = interpolate(*t,
				      betaym0, betaym1, betaym2, betaym3, betaym4);
    const double betayo = interpolate(*t,
				      betayo0, betayo1, betayo2, betayo3, betayo4);
    const double betamy = interpolate(*t,
				      betamy0, betamy1, betamy2, betamy3, betamy4);
    const double betamo = interpolate(*t,
				      betamo0, betamo1, betamo2, betamo3, betamo4);

    const double yy_infecting = (betay * Iy) / Ny;
    const double ym_infecting = (betaym * Iy) / Ny;
    const double yo_infecting = (betayo * Iy) / Ny;

    const double mm_infecting = (betam * Im) / Nm;
    const double my_infecting = (betamy * Im) / Nm;
    const double mo_infecting = (betamo * Im) / Nm;

    const double oo_infecting = (betao * Io) / No;

    const double ygot_infected = (yy_infecting + my_infecting) * Sy;
    const double ygot_infectious = a * Ey;
    const double ygot_removed = gamma * Iy;

    const double mgot_infected = (mm_infecting + ym_infecting) * Sm;
    const double mgot_infectious = a * Em;
    const double mgot_removed = gamma * Im;

    const double ogot_infected = (oo_infecting + yo_infecting + mo_infecting) * So;
    const double ogot_infectious = a * Eo;
    const double ogot_removed = gamma * Io;

    /* dSy  */ ydot[0] = -ygot_infected;
    /* dEy  */ ydot[1] = ygot_infected - ygot_infectious;
    /* dIy  */ ydot[2] = ygot_infectious - ygot_removed;
    /* dRy  */ ydot[3] = ygot_removed;

    /* dSm  */ ydot[4] = -mgot_infected;
    /* dEm  */ ydot[5] = mgot_infected - mgot_infectious;
    /* dIm  */ ydot[6] = mgot_infectious - mgot_removed;
    /* dRm  */ ydot[7] = mgot_removed;

    /* dSo  */ ydot[8] = -ogot_infected;
    /* dEo  */ ydot[9] = ogot_infected - ogot_infectious;
    /* dIo  */ ydot[10] = ogot_infectious - ogot_removed;
    /* dRo  */ ydot[11] = ogot_removed;

    yout[0] = 1/gamma * (ygot_infected + mgot_infected + ogot_infected) / (Iy + Im + Io);
    yout[1] = yout[0] * (Ny + Nm + No) / (Sy + Sm + So);

    yout[2] = 1/gamma * (yy_infecting * Sy + ym_infecting * Sm + yo_infecting * So) / Iy;
    yout[3] = 1/gamma * (my_infecting * Sy + mm_infecting * Sm + mo_infecting * So) / Im;
    yout[4] = 1/gamma *                                         (oo_infecting * So) / Io;
  }
}
