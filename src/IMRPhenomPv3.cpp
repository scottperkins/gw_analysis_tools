#include "IMRPhenomPv3.h"

/*! \file
 * IMRPhenomPv3 implementation
*/


template <class T>
void IMRPhenomPv3<T>::init_post_merger(source_parameters<T> *params, lambda_parameters<T> *lambda)
{
    this->post_merger_variables(params);
	params->f1_phase = 0.018/(params->M);
	params->f2_phase = params->fRD/2.;

	params->f1 = 0.014/(params->M);
	params->f3 = this->fpeak(params, lambda);
}

template <class T>
int IMRPhenomPv3<T>::construct_waveform(T *frequencies, /**< T array of frequencies in Hz the waveform is to be evaluated at */
				int length, /**< integer length of the array of frequencies and the waveform*/	
				std::complex<T> *waveform_plus, /** [out] complex T array for the plus polarization waveform to be output*/ 
				std::complex<T> *waveform_cross, /** [out] complex T array for the cross polarization waveform to be output */ 
				source_parameters<T> *params /* Structure of source parameters to be initialized before computation*/
				)
{
    // Imaginary unit
    const std::complex<T> I (0.,1.);
    const T half = 0.5;

    /* Store useful variables and compute derived and frequency independent variables */
    PhenomPv3Storage<T> *pv3storage;
    pv3storage = (PhenomPv3Storage<T> *)malloc(sizeof(PhenomPv3Storage<T>));
    /* Struct that stores the precession angle variables */
    sysprecquant<T> *pAngles;
    pAngles = (sysprecquant<T> *)malloc(sizeof(sysprecquant<T>));

    T deltaF = -1; // for PhenomPv3Storage
    init_PhenomPv3_Storage(pv3storage, pAngles,
        params->mass1, params->mass2,
        params->spin1x, params->spin1y, params->spin1z,
        params->spin2x, params->spin2y, params->spin2z,
        params->DL, params->incl_angle, params->phiRef,
        deltaF, frequencies[0], frequencies[length-1], params->f_ref);

    T phiHarm = 0.;
    sph_harm<T> *Y2m = (sph_harm<T> *)malloc(sizeof(sph_harm<T>));
    IMRPhenomPv3InitY2m(Y2m, params->thetaJN, phiHarm);

    // Set lambda params, should not matter to us as Pv3 is BH-only
    lambda_parameters<T> lambda;
    this->assign_lambda_param(params, &lambda);

    // Initialize post merger quantities
    this->init_post_merger(params, &lambda);
    T fCut = .2/params->M; // cut-off frequency
    int lengthCut = length; 

    useful_powers<T> pows;
	this->precalc_powers_PI(&pows);

    // Set up PN and connection coefficients 
    T deltas[6];
	T pn_amp_coeffs[7];
	T pn_phase_coeffs[12];

    this->assign_pn_amplitude_coeff(params, pn_amp_coeffs);
	this->assign_static_pn_phase_coeff(params, pn_phase_coeffs);	
	this->amp_connection_coeffs(params,&lambda,pn_amp_coeffs,deltas);
	this->phase_connection_coefficients(params,&lambda,pn_phase_coeffs);

    // Zero out the polarization arrays
    std::fill(waveform_plus, waveform_plus+length, std::complex<T> (0.));
    std::fill(waveform_cross, waveform_cross+length, std::complex<T> (0.));

    // Set up useful quantities
    const T Msec = params->M; // total mass in sec
    const T pi_Msec = GWAT_TWOPI * Msec;

    // Rescale amplitude because we aren't just using (2,2) mode anymore
	const T A0 = params->A0 * pow(Msec,7./6.) / ( 2. * sqrt(5. / (64.*M_PI)) );

    const T f_ref = params->f_ref;
	const T phic = 2*params->phi_aligned;
	const T two_pi_tc = 2*M_PI*params->tc;

    T wignerD[2][5];

    std::complex<T> PhenomDamp, expphase, half_amp_eps, hp_proj, hc_proj;
    T alpha, beta, two_epsilon, minusbeta;
    T PhenomDphase, epsphase;
    T fHz;
    for (int j = 0; j < length; j++)
    {
        fHz = frequencies[j];
        if (fHz > fCut)
        {
            lengthCut = j;
            break;
        }
        
        this->precalc_powers_ins(fHz, Msec, &pows);

        // Get non-precessing amplitude and phase
        PhenomDamp = A0 * this->build_amp(fHz, &lambda, params, &pows, pn_amp_coeffs, deltas);

        PhenomDphase = this->build_phase(fHz, &lambda, params, &pows, pn_phase_coeffs);
        PhenomDphase -= phic + two_pi_tc * (fHz - params->f_ref); // shift from initial conditions

        // Compute the precession angles
        IMRPhenomPv3_Compute_a_b_e(&alpha, &beta, &two_epsilon,
            fHz, pi_Msec, pv3storage, pAngles);
        // Compute the Wigner D elements
        minusbeta = -beta;
        IMRPhenomPv3ComputeWignerD(&wignerD, minusbeta);

        // Twist-up
        IMRPhenomPv3twist(&hp_proj, &hc_proj, alpha, Y2m, &wignerD);
        epsphase = PhenomDphase + two_epsilon;
        expphase = std::exp(-I*epsphase);
        half_amp_eps = half * PhenomDamp * expphase;

        // Store the polarizations
        waveform_plus[j] += half_amp_eps * hp_proj;
        waveform_cross[j] += half_amp_eps * hc_proj;
    }
    
    // Time correction
    T t_corr_fixed;
	if(params->shift_time)
    {
		t_corr_fixed = this->calculate_time_shift(params, &pows, pn_phase_coeffs, &lambda);
	}
	else
    {
		t_corr_fixed = 0;
	}

    // Phase correction
    T phase_corr_term;
    std::complex<T> phase_corr;
    T two_pi_t_corr = GWAT_TWOPI * t_corr_fixed;
    for (int j = 0; j<lengthCut; j++)
    {
        phase_corr_term = two_pi_t_corr * frequencies[j];
        phase_corr = std::exp(-I*phase_corr_term);

        waveform_plus[j] *= phase_corr;
        waveform_cross[j] *= phase_corr;
    }

    // Rotate by polarization angle
    std::complex<T> hp_temp, hc_temp;
    T cos_2zeta = cos(2.0 * pv3storage->zeta_polariz);
    T sin_2zeta = sin(2.0 * pv3storage->zeta_polariz);

    for (int j =0; j<lengthCut; j++)
    {
        hp_temp = waveform_plus[j];
        hc_temp = waveform_cross[j];

        waveform_plus[j] = cos_2zeta*hp_temp + sin_2zeta*hc_temp;
        waveform_cross[j] = cos_2zeta*hc_proj - sin_2zeta*hp_temp;
    }

    // Cleanup
    free(Y2m);
    free(pv3storage);
    free(pAngles);

    return 1;
}

/**
 * Compute the precession angles at a single frequency
 */
template <class T>
void IMRPhenomPv3_Compute_a_b_e(
    T *alpha, T *beta, T *two_epsilon,
    const T fHz, const T pi_Msec,
    const PhenomPv3Storage<T> *params, const sysprecquant<T> *pAngles)
{
    vector3D<T> angles;
    T xi;

    xi = pow(fHz * pi_Msec, pAngles->onethird);
    angles = compute_phiz_zeta_costhetaL3PN(xi, pAngles);

    *alpha = angles.x + params->alpha0;
    
    T epsilon = angles.y;
    *two_epsilon = 2.*epsilon;

    /* angles.z can sometimes nan just above 1.
    * The following nudge seems to be a good fix.
    */
    T one = 1;
    T diff = 1e-6;
    nudge(&(angles.z), one, diff); //This could even go from 1e-6 to 1e-15 - then would have to also change in Pv3 code. Should just fix the problem in the angles code.
    *beta = acos(angles.z);
}

/**
 * Calculate the projectors that twist-up the non-precessing mode
*/
template <class T>
void IMRPhenomPv3twist(
    std::complex<T> *hpTerm, /** [out] (twice) the plus polarization projector */
    std::complex<T> *hcTerm, /** [out] (twice) the cross polarization projector */
    const T alpha, /**< alpha angle */
    const sph_harm<T> *Y2ms, /**< spherical harmonics */
    const T (*wignerD)[2][5] /** Wigner D^2_{+/- 2, m} */
)
{
    // Imaginary unit
    const std::complex<T> I (0., 1.);
    const std::complex<T> one (1.);

    std::complex<T> Term1_sum (0.);
    std::complex<T> Term2_sum (0.);
    std::complex<T> T1, T2;

    // Computer powers of exp{i m alpha}
    std::complex<T> cexp_i_alpha = std::complex<T> (0., alpha);
    std::complex<T> cexp_2i_alpha = cexp_i_alpha*cexp_i_alpha;
    std::complex<T> cexp_mi_alpha = one/cexp_i_alpha;
    std::complex<T> cexp_m2i_alpha = cexp_mi_alpha*cexp_mi_alpha;
    std::complex<T> cexp_im_alpha[5] = {cexp_m2i_alpha, cexp_mi_alpha, one, cexp_i_alpha, cexp_2i_alpha};

    std::complex<T> Y2ms_arr[5] = 
			{Y2ms->Y2m2, Y2ms->Y2m1, Y2ms->Y20, Y2ms->Y21, Y2ms->Y22};
    std::complex<T> Y2m;

    for (int m = -2; m < 2; m++)
    {
        T1 = cexp_im_alpha[m+2] * (*wignerD)[0][m+2];
        T2 = cexp_im_alpha[-m+2] * (*wignerD)[1][m+2];

        Y2m = Y2ms_arr[m+2];
        Term1_sum += T1 * Y2m;
        Term2_sum += T2 * conj(Y2m);
    }

    *hpTerm = Term1_sum + Term2_sum;
    *hcTerm = -I*(Term1_sum - Term2_sum);
}


// Template instantiations


template void IMRPhenomPv3_Compute_a_b_e<double>(
    double *, double *, double *,
    const double, const double,
    const PhenomPv3Storage<double> *, const sysprecquant<double> *
);
template void IMRPhenomPv3_Compute_a_b_e<adouble>(
    adouble *, adouble *, adouble *,
    const adouble, const adouble,
    const PhenomPv3Storage<adouble> *, const sysprecquant<adouble> *
);

template void IMRPhenomPv3twist<double>(
    std::complex<double> *hpTerm, std::complex<double> *hcTerm,
    const double alpha, const sph_harm<double> *Y2m, const double (*wignerD)[2][5]
);
template void IMRPhenomPv3twist<adouble>(
    std::complex<adouble> *hpTerm, std::complex<adouble> *hcTerm,
    const adouble alpha, const sph_harm<adouble> *Y2m, const adouble (*wignerD)[2][5]
);

template class IMRPhenomPv3<double>;
template class IMRPhenomPv3<adouble>;
