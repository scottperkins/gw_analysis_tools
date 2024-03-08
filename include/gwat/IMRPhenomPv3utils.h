#ifndef IMRPHENOMPV3_UTILS_H
#define IMRPHENOMPV3_UTILS_H

#include "util.h"
#include <stdexcept>
#include <gsl/gsl_sf.h>
#include <vector>
#include "IMRPhenomP.h"


typedef std::vector<double> DOUBLEVector;


/*! \struct
 * 3D vector
 */
typedef struct tagvector
{
    double x;
    double y;
    double z;
} vector;

/*! \struct
 * Holds quantities needed for precession calculations.
 */
typedef struct tagsystemprecessionquantities
{
    double onethird;
    double constants_u[6];
    double constants_phiz[6];
    double constants_zeta[6];
    double constants_L[6];
    double phiz_0, zeta_0, constant_of_S;
    double c_1, Ssqave, sqrtSsqave, Seff, c1_2, nu_2, nu_4, c_1_over_nu, S0_norm;
    double S1_norm_2, S2_norm_2;
    double dot1, dot2, dot12, dot1n, dot2n;
    double deltam_over_M, nu, q;
} sysprecquant;

/**
 * Structure storing initial and derived variables for IMRPhenomPv3
 */
typedef struct tagPhenomPv3Storage
{
    int PRECESSING;   /**< integer to signify if system is precessing, 1 for false (not precessing), 0 for true (precessing) */
    /* input parameters */
    double m1_SI;       /**< mass of primary in SI (kg) */
    double m2_SI;       /**< mass of secondary in SI (kg) */
    double chi1x;       /**< x-component of dimensionless spin on primary w.r.t. Lhat = (0,0,1) */
    double chi1y;       /**< y-component of dimensionless spin on primary w.r.t. Lhat = (0,0,1) */
    double chi1z;       /**< z-component of dimensionless spin on primary w.r.t. Lhat = (0,0,1) */
    double chi2x;       /**< x-component of dimensionless spin on secondary w.r.t. Lhat = (0,0,1) */
    double chi2y;       /**< y-component of dimensionless spin on secondary w.r.t. Lhat = (0,0,1) */
    double chi2z;       /**< z-component of dimensionless spin on secondary w.r.t. Lhat = (0,0,1) */
    double distance_SI; /**< distance to source in SI (m) */
    double inclination; /**< inclination - used to compute the angle thetaJN (rad) */
    double phiRef;      /**< */
    double deltaF;      /**< frequency spacing (Hz) */
    double f_min;       /**< starting GW frequency (Hz) */
    double f_max;       /**< ending GW frequency (Hz) */
    double f_ref;       /**< reference GW frequency (Hz) */
    /* derived parameters */
    double m1_Msun;      /**< mass of primary in solar masses */
    double m2_Msun;      /**< mass of secondary in solar masses */
    double Mtot_SI;      /**< total mass in SI (kg) */
    double Mtot_Msun;    /**< total mass in solar masses */
    double eta;          /**< Symmetric mass ratio*/
    double q;            /* with m1>=m2 so q>=1 */
    double Msec;         /**< Total mass in seconds */
    double f_ref_Orb_Hz; /**< Reference orbital frequency (Hz) [It's the reference GW frequency converted to orbital frequency] */
    double twopi_Msec;   /**< LAL_TWOPI * Msec */
    double amp0;         /**< frequency domain physical scaling */
    /* variables used when rotating input parameters (LAL frame) into PhenomP intrinsic parameters  */
    double chip;         /**< effective precessing parameter */
    double thetaJN;      /**< Angle between J0 and line of sight (z-direction) */
    double alpha0;       /**< Initial value of alpha angle (azimuthal precession angle) */
    double phi_aligned;  /**< Initial phase to feed the underlying aligned-spin model */
    double zeta_polariz; /**< Angle to rotate the polarizations */
    /* compute spins in polar coordinates */
    double chi1_mag;   /**< dimensionless spin magnitude on primary */
    double chi1_theta; /**< polar angle w.r.t. Lhat = (0,0,1) on primary */
    double chi1_phi;   /**< azimuthal angle w.r.t. Lhat = (0,0,1) on primary */
    double chi2_mag;   /**< dimensionless spin magnitude on secondary */
    double chi2_theta; /**< polar angle w.r.t. Lhat = (0,0,1) on secondary */
    double chi2_phi;   /**< azimuthal angle w.r.t. Lhat = (0,0,1) on secondary */
    /* Precession angles at reference frequency */
    double alphaRef;   /**< azimuthal precession angle at f_ref */
    double epsilonRef; /**< epsilon precession angle at f_ref */
    double betaRef;    /**< beta (opening angle) precession angle at f_ref */
} PhenomPv3Storage;

void PhenomPrecessingSpinEnforcePrimary(double *m1, double *m2,
    double *chi1x, double *chi1y, double *chi1z,
    double *chi2x, double *chi2y,double *chi2z
);

void nudge(double *x, double X, double epsilon);

static vector CreateSphVector(const double r, const double th, const double ph);
static vector ScaleVector(double c, vector vec);
static vector VectorSum(vector vec1, vector vec2);
static double DotProduct(const vector vec1, const vector vec2);
static double VectorNorm(const vector vec);
static vector CrossProduct(const vector vec1, const vector vec2);

static vector Roots(const double L_norm, const double J_norm, const sysprecquant *system);
static vector BCDcoeff(const double L_norm, const double J_norm, const sysprecquant *system);
static double beta(const double a, const double b, const sysprecquant *system);
static double sigma(const double a, const double b, const sysprecquant *system);
static double tau(const double a, const double b, const sysprecquant *system);

static double L_norm_3PN_of_xi(const double xi, const double xi_2, const double L_norm, const sysprecquant *system);
static double J_norm_of_xi(const double L_norm, const sysprecquant *system);
static double S_norm_of_xi(const double xi, const double xi_2, const vector roots, const sysprecquant *system);
static double costhetaL(const double J_norm, const double L_norm, const double S_norm);

static double u_of_xi(const double xi, const double xi_2, const sysprecquant *system);
static double phiz_of_xi(const double xi, const double xi_2, const double J_norm, const sysprecquant *system);
static double zeta_of_xi(const double xi, const double xi_2, const sysprecquant *system);

static vector computeMScorrections(const double xi, const double xi_2, const double L_norm, const double J_norm, const vector roots, const sysprecquant *system);

static int checkOmegaz5(const double Omegaz5);
static void invalidExpansionOrder(const int ExpansionOrder);

static vector compute_phiz_zeta_costhetaL3PN(const double xi, const sysprecquant *system);

static double L2PNR(const double v, const double eta);
double L3PN(
    const double f_orb_hz,
    const double m1, const double m2,
    const double s1x, const double s1y, const double s1z,
    const double s2x, const double s2y, const double s2z,
    const double f_0, const int ExpansionOrder);

void CartesianToPolar(double *polar, double *azimuthal, double *magnitude, double x, double y, double z);

void PhenomP_Param_Transform(
    double *chi1_l, double *chi2_l, double *chip, double *thetaJN, double *alpha0, double *phi_aligned, double *zeta_polariz,
    const double m1, const double m2, const double f_ref, const double phiRef, const double incl,
    const double s1x, const double s1y, const double s1z, const double s2x, const double s2y, const double s2z,
    IMRPhenomP_version_type IMRPhenomP_version);

static void InitializePrecession(
    sysprecquant* system,
    const double m1, const double m2, const double mul, const double phl, const double mu1, const double ph1,
    const double ch1, const double mu2, const double ph2, const double ch2, const double f_0,
    const int ExpansionOrder);

void OrbitalAngMom3PNSpinning(
    DOUBLEVector *L_norm_3PN, DOUBLEVector *f_orb_hz,
    const double m1_SI, const double m2_SI,
    const double mul, const double phl, double mu1, double ph1, double ch1, double mu2, double ph2, double ch2,
    const double f_0, const int ExpansionOrder);


#endif /* IMRPHENOMPV3_UTILS_H */