#include "IMRPhenomPv3.h"

/*! \file
 * IMRPhenomPv3 implementation
*/


/**
 * Precomputes useful quantities and populates the
 * PhenomPv3HMStorage and sysq (for precession angles) structs.
 */
template <class T> static int init_PhenomPv3_Storage(
    PhenomPv3Storage<T> *p,   /**< [out] PhenomPv3Storage struct */
    sysprecquant<T> *pAngles,           /**< [out] precession angle pre-computations struct */
    T m1_SI,             /**< mass of primary in SI (kg) */
    T m2_SI,             /**< mass of secondary in SI (kg) */
    T S1x,               /**< x-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    T S1y,               /**< y-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    T S1z,               /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    T S2x,               /**< x-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    T S2y,               /**< y-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    T S2z,               /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    const T distance,    /**< distance of source (m) */
    const T inclination, /**< inclination of source (rad) */
    const T phiRef,      /**< reference orbital phase (rad) */
    const T deltaF,      /**< Sampling frequency (Hz) */
    const T f_min,       /**< Starting GW frequency (Hz) */
    const T f_max,       /**< End frequency; 0 defaults to ringdown cutoff freq */
    const T f_ref        /**< Reference GW frequency (Hz) */
)
{
    p->PRECESSING = 0;
    if (S1x == 0. && S1y == 0. && S2x == 0. && S2y == 0.)
    {
        p->PRECESSING = 1; // This means the system is not precessing
    }

    /* input parameters */
    p->m1_SI = m1_SI;
    p->m2_SI = m2_SI;
    p->chi1x = S1x;
    p->chi1y = S1y;
    p->chi1z = S1z;
    p->chi2x = S2x;
    p->chi2y = S2y;
    p->chi2z = S2z;
    p->distance_SI = distance;
    p->phiRef = phiRef;
    p->deltaF = deltaF;
    p->f_min = f_min;
    p->f_max = f_max;
    p->f_ref = f_ref;

    PhenomPrecessingSpinEnforcePrimary(
        &(p->m1_SI),
        &(p->m2_SI),
        &(p->chi1x),
        &(p->chi1y),
        &(p->chi1z),
        &(p->chi2x),
        &(p->chi2y),
        &(p->chi2z));

    p->m1_Msun = m1_SI / GWAT_MSUN_SI;
    p->m2_Msun = m2_SI / GWAT_MSUN_SI;
    p->Mtot_SI = p->m1_SI + p->m2_SI;
    p->Mtot_Msun = p->m1_Msun + p->m2_Msun;

    p->eta = p->m1_Msun * p->m2_Msun / (p->Mtot_Msun * p->Mtot_Msun);
    p->q = p->m1_Msun / p->m2_Msun; /* with m1>=m2 so q>=1 */

    /* check for rounding errors */
    if (p->eta > 0.25 || p->q < 1.0)
    {
        nudge(&(p->eta), 0.25, 1e-6);
        nudge(&(p->q), 1.0, 1e-6);
    }

    p->Msec = p->Mtot_Msun * GWAT_MTSUN_SI; /* Total mass in seconds */

    p->amp0 = (p->Mtot_Msun) * GWAT_MRSUN_SI * (p->Mtot_Msun) * GWAT_MTSUN_SI / (p->distance_SI);

    /* Rotate to PhenomP frame */
    /* chi1_l == chi1z, chi2_l == chi2z for intermediate calculations*/
    T chi1_l, chi2_l;
    PhenomP_ParametersFromSourceFrame(
        &chi1_l, &chi2_l, &(p->chip), &(p->thetaJN), &(p->alpha0), &(p->phi_aligned), &(p->zeta_polariz),
        p->m1_SI, p->m2_SI, p->f_ref, p->phiRef, inclination,
        p->chi1x, p->chi1y, p->chi1z,
        p->chi2x, p->chi2y, p->chi2z, IMRPhenomPv3_V);

    p->inclination = p->thetaJN;

    if (p->PRECESSING != 1) // precessing case. compute angles
    {
        /* Initialize precession angles */
        /* evaluating the angles at the reference frequency */
        p->f_ref_Orb_Hz = 0.5 * p->f_ref; /* factor of 0.5 to go from GW to Orbital frequency */

        /* precompute everything needed to compute precession angles */
        /* ExpansionOrder specifies how many terms in the PN expansion of the precession angles to use.
        * In PhenomP3 we set this to 5, i.e. all but the highest order terms.
        * */
        int ExpansionOrder = 5;
        InitializePrecession(
            pAngles,
            p->m1_SI, p->m2_SI,
            1.0, 0.0,
            cos(p->chi1_theta), p->chi1_phi, p->chi1_mag,
            cos(p->chi2_theta), p->chi2_phi, p->chi2_mag,
            p->f_ref, ExpansionOrder);
    }
}

/**
 * Compute the precession angles at a single frequency
 */
template <class T>
static void IMRPhenomPv3_Compute_a_b_e(
    T *alpha, T *beta, T *two_epsilon,
    T fHz, const T pi_Msec,
    PhenomPv3Storage<T> *params, sysprecquant<T> *pAngles)
{
    vector3D<T> angles;
    T xi;

    xi = pow(fHz * pi_Msec, pAngles->onethird);
    angles = compute_phiz_zeta_costhetaL3PN(xi, pAngles);

    *alpha = angles.x + params->alpha0;
    
    T epsilon = angles.y
    *two_epsilon = 2.*epsilon;

    /* angles.z can sometimes nan just above 1.
    * The following nudge seems to be a good fix.
    */
    nudge(&(angles.z), 1.0, 1e-6); //This could even go from 1e-6 to 1e-15 - then would have to also change in Pv3 code. Should just fix the problem in the angles code.
    *beta = acos(angles.z);
}
