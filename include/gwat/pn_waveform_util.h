#ifndef PN_WAVEFORM_UTIL_H
#define PN_WAVEFORM_UTIL_H
/*! \file
 * Header file for all PN-waveform-specific utilities
 */
template<class T>
T t_2PN(T f, T eta, T chirpmass, T chi1, T chi2, T tc);

template<class T>
T t_0PN(T f, T chirpmass);
template<class T>
T f_0PN(T t, T chirpmass);

template<class T>
T FISCO(T mass);
#endif
