#ifndef IMRPHENOMD_NRT_EA_H
#define IMRPHENOMD_NRT_EA_H
#include "IMRPhenomD.h"
#include "IMRPhenomD_NRT.h"
#include "util.h"
#include "waveform_generator.h"

template<class T>
class EA_IMRPhenomD_NRT: public IMRPhenomD_NRT<T>
{
public:
  virtual T EA_phase_ins(T f, useful_powers<T> *powers, source_parameters<T> *param);

  virtual T EA_amp_ins(T f, useful_powers<T> *powers, source_parameters<T> *param);

  virtual int EA_construct_waveform(T *frequencies, int length, waveform_polarizations<T> *waveform, source_parameters<T> *params); 

  virtual int construct_phase(T *frequencies, int length, T *phase, source_parameters<T> *params);
  
  virtual int construct_amplitude(T *frequencies, int length, T *amplitude, source_parameters<T> *params);
  
  virtual int construct_waveform(T *frequencies, int length, std::complex<T> *waveform, source_parameters<T> *params);
};

#endif
