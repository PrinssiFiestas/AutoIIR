/* MIT License
 * Copyright (c) 2023 Lauri Fiestas
 * https://github.com/PrinssiFiestas/AutoIIR/blob/main/LICENSE */

#ifndef AUTOIIR_H
#define AUTOIIR_H

#ifdef __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
//
//			COMPILE-TIME OPTIONS
//
//-----------------------------------------------------------------------------
// Set these flags globally
// autoMem filters use static memory to store filter state so the amount of
// of filter calls, audio inputs and poles are finite. Increase when running
// out and decreace if limited memory is an issue.

#ifndef AUTOIIR_MAX_FILTERS
#define AUTOIIR_MAX_FILTERS 200
#endif

#ifndef AUTOIIR_MAX_AUDIO_INPUTS
#define AUTOIIR_MAX_AUDIO_INPUTS 24 // 22 + 2 = max channels supported by VST3
#endif

#ifndef AUTOIIR_MAX_POLES
#define AUTOIIR_MAX_POLES 16
#endif

// By default there is no error handling to save performance. This may cause
// the program to crash on errors. By setting this flag error messages are sent
// to stderr and filter functions return without crashing. User can also get
// access to error messages with callback set by iir_setErrorMessageCallback().
//#define AUTOIIR_DEBUG

// Set this flag for files with namespacing problems
//#define AUTOIIR_NO_SHORT_ALIASES

//-----------------------------------------------------------------------------
//
//			PUBLIC API
//
//-----------------------------------------------------------------------------

#include <stddef.h>

#ifndef PI
#define PI 3.14159265359
#endif

#ifndef SQRT2
#define SQRT2 1.41421356237
#endif

#ifndef INV_SQRT2
#define INV_SQRT2 0.707106781187
#endif

// Set internal sample rate for filters that can take frequency in hertz
// Optional if normalized frequencies from 0.0 to 1.0 are used.
void iir_setSampleRate(double);

// Resets filter counting and sets channel for autoMem filters
// Should be the first function called in processing loop.
// Optional if no autoMem filters are used.
// Channel can just be set to 0 when processing in mono.
void autoMem_initProcessing(unsigned int channel);

// Filter state
// Use for manual memory management for manual_filter().
typedef struct iir_Memory
{
	double x[AUTOIIR_MAX_POLES + 1];
	double y[AUTOIIR_MAX_POLES + 1];
} iir_Memory;

// Coefficients for the recursion equation
// This library calculates filter output with recursion equation in form of
// 		y0 = a0x0 + a1x1 + ... + anxn + b1y1 + b2y2 + ... + bnyn
// 	as opposed to difference equation in form of
// 		y0 = a0x0 + a1x1 + ... + anxn - b1y1 - b2y2 - ... - bnyn
// 	so note the signs of b coefficients when writing your own filters. 
typedef struct iir_Coeffs
{
	double a[AUTOIIR_MAX_POLES + 1];
	double b[AUTOIIR_MAX_POLES + 1];
	unsigned int poles;
} iir_Coeffs;

// Table of non-linear functions for non-linear filters
typedef struct iir_FuncTable
{
	double(*fx[AUTOIIR_MAX_POLES + 1])(double);
	double(*fy[AUTOIIR_MAX_POLES + 1])(double);
} iir_FuncTable;

// Returns table with all functions set to a functions that do nothing
iir_FuncTable initFuncTable(void);

// Returns table with all functions set to function passed as argument
iir_FuncTable initFuncTableWith(double(*func)(double));

#define IIR_LINEAR NULL

double autoMem_filter(double input, iir_Coeffs*, iir_FuncTable* nonLinearities);
double manual_filter(iir_Memory*, double input, iir_Coeffs*, iir_FuncTable* nonLinearities);

// Sets all autoMem filter states to 0.0
void autoMem_clear(void);

// Sets filter state to 0.0
void Memory_clear(iir_Memory* me);

// Set callback to be called on errors when AUTOIIR_DEBUG is defined
// By default error messages are sent to stderr. If AUTOIIR_DEBUG isn't defined
// this function does nothing and ther's no error handling to save performance.
void iir_setErrorMessageCallback(void(*callback)(const char* errorMessage));

//-----------------------------------------------------------------------------

#ifndef AUTOIIR_NO_SHORT_ALIASES

// Accurate filters
// Frequencies between 0.0 to 1.0 are assumed to be normalized.
// Frequencies greater than 1.0 are assumed to be in hertz.
// Set sample rate with iir_setSampleRate() when using hertz!

#define lowPass6(input, freq)               autoMem_lowPass6(input, freq)
#define highPass6(input, freq)              autoMem_highPass6(input, freq)
#define lowPass12(input, freq, q)           autoMem_lowPass12(input, freq, q)
#define highPass12(input, freq, q)          autoMem_highPass12(input, freq, q)
#define bandPass12(input, freq, q)          autoMem_bandPass12(input, freq, q)
#define bandStop12(input, freq, q)          autoMem_bandStop12(input, freq, q)
#define butterworthLp12(input, freq)        autoMem_butterworthLp12(input, freq)
#define butterworthHp12(input, freq)        autoMem_butterworthHp12(input, freq)
#define allPass6(input, freq)               autoMem_allPass6(input, freq)
#define allPass12(input, freq, q)           autoMem_allPass12(input, freq, q)
#define lowShelving(input, freq, gain_dB)   autoMem_lowShelving(input, freq, gain_dB)
#define highShelving(input, freq, gain_dB)  autoMem_highShelving(input, freq, gain_dB)
#define peak(input, freq, q, gain_dB)       autoMem_peak(input, freq, q, gain_dB)
#define peakConstQ(input, freq, q, gain_dB) autoMem_peakConstQ(input, freq, q, gain_dB)
#define linkwitzRileyLp12(input, freq)      autoMem_linkwitzRileyLp12(input, freq)
#define linkwitzRileyHp12(input, freq)      autoMem_linkwitzRileyHp12(input, freq)

// Speed optimized filters
// Only takes normalized frequencies from 0.0 to 1.0 excluding 0.0 and 1.0.
// Frequency and Q won't be accurate.

#define fastLp6(input, freq)                    autoMem_fastLp6(input, freq)
#define fastHp6(input, freq)                    autoMem_fastHp6(input, freq)
#define fastLp12(input, freq, damp)             autoMem_fastLp12(input, freq, damp)
#define fastHp12(input, freq, damp)             autoMem_fastHp12(input, freq, damp)
// TODO
//#define fastBp12(input, freq, q)                autoMem_fastBp12(input, freq, q)
//#define fastBs12(input, freq, q)                autoMem_fastBs12(input, freq, q)
//#define fastButtLp12(input, freq)               autoMem_fastButtLp12(input, freq)
//#define fastButtHp12(input, freq)               autoMem_fastButtHp12(input, freq)
//#define fastAp6(input, freq)                    autoMem_fastAp6(input, freq)
//#define fastAp12(input, freq, bandWidth)        autoMem_fastAp12(input, freq, bandWidth)
//#define fastLShelving(input, freq)              autoMem_fastLShelving(input, freq)
//#define fastHShelving(input, freq)              autoMem_fastHShelving(input, freq)
//#define fastPeak(input, freq, q, gain)          autoMem_fastPeak(input, freq, q, gain)
//#define fastPeakConstQ(input, freq, q, gain)    autoMem_fastPeakConstQ(input, freq, q, gain)
//#define fastLinkwitzRileyLp12(input, freq)      autoMem_fastLinkwitzRileyLp12(input, freq)
//#define fastLinkwitzRileyHp12(input, freq)      autoMem_fastLinkwitzRileyHp12(input, freq)

#endif // AUTOIIR_NO_SHORT_ALIASES

//-----------------------------------------------------------------------------

double autoMem_lowPass6(double input, double freq);
double autoMem_highPass6(double input, double freq);
double autoMem_lowPass12(double input, double freq, double q);
double autoMem_highPass12(double input, double freq, double q);
double autoMem_bandPass12(double input, double freq, double q);
double autoMem_bandStop12(double input, double freq, double q);
double autoMem_butterworthLp12(double input, double freq);
double autoMem_butterworthHp12(double input, double freq);
double autoMem_allPass6(double input, double freq);
double autoMem_allPass12(double input, double freq, double q);
double autoMem_lowShelving(double input, double freq, double gain_dB);
double autoMem_highShelving(double input, double freq, double gain_dB);
double autoMem_peak(double input, double freq, double q, double gain_dB);
double autoMem_peakConstQ(double input, double freq, double q, double gain_dB);
double autoMem_linkwitzRileyLp12(double input, double freq);
double autoMem_linkwitzRileyHp12(double input, double freq);

double autoMem_fastLp6(double input, double normalizedFreq);
double autoMem_fastHp6(double input, double normalizedFreq);
double autoMem_fastLp12(double input, double normalizedFreq, double damping);
double autoMem_fastHp12(double input, double normalizedFreq, double damping);
// TODO
//double autoMem_fastBp12(double input, double freq, double q);
//double autoMem_fastBs12(double input, double freq, double q);
//double autoMem_fastButtLp12(double input, double freq);
//double autoMem_fastButtHp12(double input, double freq);
//double autoMem_fastAp6(double input, double freq);
//double autoMem_fastAp12(double input, double freq, double bandWidth);
//double autoMem_fastLShelving(double input, double freq);
//double autoMem_fastHShelving(double input, double freq);
//double autoMem_fastPeak(double input, double freq, double q, double gain);
//double autoMem_fastPeakConstQ(double input, double freq, double q, double gain);
//double autoMem_fastLinkwitzRileyLp12(double input, double freq);
//double autoMem_fastLinkwitzRileyHp12(double input, double freq);

//-----------------------------------------------------------------------------

iir_Coeffs Coeffs_lowPass6(double freq);
iir_Coeffs Coeffs_highPass6(double freq);
iir_Coeffs Coeffs_lowPass12(double freq, double q);
iir_Coeffs Coeffs_highPass12(double freq, double q);
iir_Coeffs Coeffs_bandPass12(double freq, double q);
iir_Coeffs Coeffs_bandStop12(double freq, double q);
iir_Coeffs Coeffs_butterworthLp12(double freq);
iir_Coeffs Coeffs_butterworthHp12(double freq);
iir_Coeffs Coeffs_allPass6(double freq);
iir_Coeffs Coeffs_allPass12(double freq, double bandWidth);
iir_Coeffs Coeffs_lowShelving(double freq, double gain);
iir_Coeffs Coeffs_highShelving(double freq, double gain_dB);
iir_Coeffs Coeffs_peak(double freq, double q, double gain);
iir_Coeffs Coeffs_peakConstQ(double freq, double q, double gain);
iir_Coeffs Coeffs_linkwitzRileyLp12(double freq);
iir_Coeffs Coeffs_linkwitzRileyHp12(double freq);

iir_Coeffs Coeffs_fastLp6(double normalizedFreq);
iir_Coeffs Coeffs_fastHp6(double normalizedFreq);
iir_Coeffs Coeffs_fastLp12(double normalizedFreq, double dampingFactor);
iir_Coeffs Coeffs_fastHp12(double normalizedFreq, double dampingFactor);
// TODO
//iir_Coeffs Coeffs_fastBp12(double freq, double q);
//iir_Coeffs Coeffs_fastBs12(double freq, double q);
//iir_Coeffs Coeffs_fastButtLp12(double freq);
//iir_Coeffs Coeffs_fastButtHp12(double freq);
//iir_Coeffs Coeffs_fastAp6(double freq);
//iir_Coeffs Coeffs_fastAp12(double freq, double bandWidth);
//iir_Coeffs Coeffs_fastLShelving(double freq);
//iir_Coeffs Coeffs_fastHShelving(double freq);
//iir_Coeffs Coeffs_fastPeak(double freq, double q, double gain);
//iir_Coeffs Coeffs_fastPeakConstQ(double freq, double q, double gain);
//iir_Coeffs Coeffs_fastLinkwitzRileyLp12(double freq);
//iir_Coeffs Coeffs_fastLinkwitzRileyHp12(double freq);

//-----------------------------------------------------------------------------
//
//			IMPLEMENTATION
//
//-----------------------------------------------------------------------------

#ifdef AUTOIIR_IMPLEMENTATION

#include <math.h>

#ifdef AUTOIIR_DEBUG
#include <stdio.h>
static void(*iir_errorMessageCallback)(const char*) = perror;
#endif

static unsigned int autoMem_filterCallCount;
static unsigned int autoMem_channel;
static iir_Memory autoMem[AUTOIIR_MAX_FILTERS][AUTOIIR_MAX_AUDIO_INPUTS];
static double iir_sampleRate = 1./*to prevent zero division*/;

void iir_setErrorMessageCallback(void(*callback)(const char* message))
{
	#ifdef AUTOIIR_DEBUG
	iir_errorMessageCallback = callback;
	#else // get rid of unused variable warning
	(void)callback;
	#endif
}

void iir_setSampleRate(double sampleRate)
{
	#ifdef AUTOIIR_DEBUG
	if (sampleRate <= 0.)
	{
		iir_errorMessageCallback(
				"Sample rate passed to iir_setSampleRate() must be greater than 0!");
		return;
	}
	#endif
	iir_sampleRate = sampleRate;
}

void autoMem_initProcessing(unsigned int channel)
{
	autoMem_filterCallCount = 0;
	#ifdef AUTOIIR_DEBUG
	if (channel > AUTOIIR_MAX_AUDIO_INPUTS)
	{
		iir_errorMessageCallback("Exceeded AUTOIIR_MAX_AUDIO_INPUTS!");
		channel = 0;
		return;
	}
	#endif
	autoMem_channel = channel;
}

double autoMem_filter(double input, iir_Coeffs* coeffs, iir_FuncTable* nonLinearities)
{
	iir_Memory* uniqueMemory = &autoMem[autoMem_filterCallCount][autoMem_channel];
	autoMem_filterCallCount++;
	#ifdef AUTOIIR_DEBUG
	if (autoMem_filterCallCount > AUTOIIR_MAX_FILTERS)
	{
		iir_errorMessageCallback("Exceeded AUTOIIR_MAX_FILTERS!");
		autoMem_filterCallCount = 0;
		return 0.;
	}
	#endif
	return manual_filter(uniqueMemory, input, coeffs, nonLinearities);
}

double manual_filter(iir_Memory* mem, double x, iir_Coeffs* coeff, iir_FuncTable* nonLins)
{
	mem->x[0] = x;
	mem->y[0] = coeff->a[0]*x;
	#ifdef AUTOIIR_DEBUG
	if (coeff->poles > AUTOIIR_MAX_POLES)
	{
		iir_errorMessageCallback("Exceeded AUTOIIR_MAX_POLES!");
		return 0.;
	}
	#endif
	const unsigned int POLES = coeff->poles == 0 ? AUTOIIR_MAX_POLES : coeff->poles;
	for(unsigned int i = POLES; i > 0; i--)
	{
		// Calculate recursion equation
		// Note signs when writing coefficients!
		if (nonLins)
			mem->y[0] += coeff->a[i]*nonLins->fx[i](mem->x[i])
		    	       + coeff->b[i]*nonLins->fy[i](mem->y[i]);
		else
			mem->y[0] += coeff->a[i]*mem->x[i] + coeff->b[i]*mem->y[i];

		// Update filter memory
		mem->x[i] = mem->x[i - 1];
		mem->y[i] = mem->y[i - 1];
	}

	const double antiDenormal = 1e-30;
	mem->y[0] += antiDenormal;
	mem->y[0] -= antiDenormal;

	return mem->y[0];
}

static double returnArgWithoutDoingAnything(double arg)
{
	return arg;
}

iir_FuncTable initFuncTable(void)
{
	return initFuncTableWith(returnArgWithoutDoingAnything);
}

iir_FuncTable initFuncTableWith(double(*func)(double))
{
	iir_FuncTable table;
	for (int i = 0; i < AUTOIIR_MAX_POLES; i++)
		table.fx[i] = table.fy[i] = func;
	return table;
}

void autoMem_clear(void)
{
	for (int i = 0; i < AUTOIIR_MAX_FILTERS; i++)
		for (int j = 0; j < AUTOIIR_MAX_AUDIO_INPUTS; j++)
			Memory_clear(&autoMem[i][j]);
}

void Memory_clear(iir_Memory* me)
{
	for (int i = 0; i < AUTOIIR_MAX_POLES; i++)
		me->x[i] = me->y[i] = 0.;
}

//-----------------------------------------------------------------------------
//
//			FILTERS
//
//-----------------------------------------------------------------------------
//
// TODO
// Clamping frequencies for all non-fast filters
//
//-----------------------------------------------------------------------------
//		LOW PASS 6

iir_Coeffs Coeffs_lowPass6(double freq)
{
	const unsigned int poles = 1;

	double freq1 = freq <= 1. ? freq : 2.*freq/iir_sampleRate;
	double gamma = 2. - cos(PI*freq1);
	double b1 = gamma - sqrt(gamma*gamma - 1.);
	double a0 = 1. - b1;

	iir_Coeffs coeff = {{a0},{0., b1}, poles};
	return coeff;
}

double autoMem_lowPass6(double input, double freq)
{
	iir_Coeffs coeffs = Coeffs_lowPass6(freq);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		HIGH PASS 6

iir_Coeffs Coeffs_highPass6(double freq)
{
	const unsigned int poles = 1;

	double freq1 = freq <= 1. ? freq : 2.*freq/iir_sampleRate;
	double theta = PI*freq1;
	double gamma = cos(theta)/(1. + sin(theta));
	double a0 = (1. + gamma)/2.;
	double a1 = -a0;
	double b1 = gamma;

	iir_Coeffs coeff = {{a0, a1},{0., b1}, poles};
	return coeff;
}

double autoMem_highPass6(double input, double freq)
{
	iir_Coeffs coeffs = Coeffs_highPass6(freq);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		LOW PASS 12

iir_Coeffs Coeffs_lowPass12(double freq, double q)
{
	const unsigned int poles = 2;

	double freq1 = freq <= 1. ? freq : 2.*freq/iir_sampleRate;
	double theta = PI*freq1;
	double d     = 1./q;
	double dsint = d*sin(theta)/2.;
	double beta  = .5*(1. - dsint)/(1. + dsint);
	double gamma = (.5 + beta)*cos(theta);
	double delta = .5 + beta - gamma;
	double a0 = delta/2.;
	double a1 = delta;
	double a2 = a0;
	double b1 = 2.*gamma;
	double b2 = -2.*beta;
	
	iir_Coeffs coeff = {{a0, a1, a2},{0., b1, b2}, poles};
	return coeff;
}

double autoMem_lowPass12(double input, double freq, double q)
{
	iir_Coeffs coeffs = Coeffs_lowPass12(freq, q);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		HIGH PASS 12

iir_Coeffs Coeffs_highPass12(double freq, double q)
{
	const unsigned int poles = 2;

	double freq1 = freq <= 1. ? freq : 2.*freq/iir_sampleRate;
	double theta = PI*freq1;
	double d     = 1./q;
	double dsint = d*sin(theta)/2.;
	double beta  = .5*(1. - dsint)/(1. + dsint);
	double gamma = (.5 + beta)*cos(theta);
	double sigma = .5 + beta + gamma;
	double a0 = sigma/2.;
	double a1 = -sigma;
	double a2 = a0;
	double b1 = 2.*gamma;
	double b2 = -2.*beta;
	
	iir_Coeffs coeff = {{a0, a1, a2},{0., b1, b2}, poles};
	return coeff;
}

double autoMem_highPass12(double input, double freq, double q)
{
	iir_Coeffs coeffs = Coeffs_highPass12(freq, q);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		BAND PASS 12

iir_Coeffs Coeffs_bandPass12(double freq, double q)
{
	const unsigned int poles = 2;

	double freq1 = freq <= 1. ? freq/2. : freq/iir_sampleRate;
	double k = tan(PI*freq1);
	double delta = k*k*q + k + q;
	double a0 = k/delta;
	double a1 = 0.;
	double a2 = -a0;
	double b1 = -2.*q*(k*k - 1.)/delta;
	double b2 = -(k*k*q - k + q)/delta;
	
	iir_Coeffs coeff = {{a0, a1, a2},{0., b1, b2}, poles};
	return coeff;
}

double autoMem_bandPass12(double input, double freq, double q)
{
	iir_Coeffs coeffs = Coeffs_bandPass12(freq, q);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		BAND STOP 12

iir_Coeffs Coeffs_bandStop12(double freq, double q)
{
	const unsigned int poles = 2;

	double freq1 = freq <= 1. ? freq/2. : freq/iir_sampleRate;
	double k = tan(PI*freq1);
	double delta = k*k*q + k + q;
	double a0 = q*(k*k + 1.)/delta;
	double a1 = 2.*q*(k*k - 1.)/delta;
	double a2 = a0;
	double b1 = -a1;
	double b2 = -(k*k*q - k + q)/delta;
	
	iir_Coeffs coeff = {{a0, a1, a2},{0., b1, b2}, poles};
	return coeff;
}

double autoMem_bandStop12(double input, double freq, double q)
{
	iir_Coeffs coeffs = Coeffs_bandStop12(freq, q);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		BUTTERWORTH LOW PASS 12

iir_Coeffs Coeffs_butterworthLp12(double freq)
{
	const unsigned int poles = 2;

	double freq1 = freq <= 1. ? freq/2. : freq/iir_sampleRate;
	double bounded = fmin(fmax(freq1, 1e-5), .5 - 1e-4); // prevent singularities
	double c = 1./tan(PI*bounded);
	double a0 = 1./(1. + SQRT2*c + c*c);
	double a1 = 2.*a0;
	double a2 = a0;
	double b1 = -2.*a0*(1. - c*c);
	double b2 = -a0*(1. - SQRT2*c + c*c);
	
	iir_Coeffs coeff = {{a0, a1, a2},{0., b1, b2}, poles};
	return coeff;
}

double autoMem_butterworthLp12(double input, double freq)
{
	iir_Coeffs coeffs = Coeffs_butterworthLp12(freq);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		BUTTERWORTH HIGH PASS 12

iir_Coeffs Coeffs_butterworthHp12(double freq)
{
	const unsigned int poles = 2;

	double freq1 = freq <= 1. ? freq/2. : freq/iir_sampleRate;
	double _freq1 = .5 - freq1;
	double bounded = fmin(fmax(_freq1, 1e-5), .5 - 1e-4); // prevent singularities
	double c = 1./tan(PI*bounded);
	double a0 = 1./(1. + SQRT2*c + c*c);
	double a1 = -2.*a0;
	double a2 = a0;
	double b1 = -2.*a0*(c*c - 1.);
	double b2 = -a0*(1. - SQRT2*c + c*c);

	iir_Coeffs coeff = {{a0, a1, a2},{0., b1, b2}, poles};
	return coeff;
}

double autoMem_butterworthHp12(double input, double freq)
{
	iir_Coeffs coeffs = Coeffs_butterworthHp12(freq);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		ALL PASS 6

iir_Coeffs Coeffs_allPass6(double freq)
{
	const unsigned int poles = 1;

	double freq1 = freq <= 1. ? freq/2. : freq/iir_sampleRate;
	double theta = PI*freq1;
	double alpha = (tan(theta) - 1.)/(tan(theta) + 1.);
	double a0 = alpha;
	double a1 = 1.;
	double b1 = -alpha;

	iir_Coeffs coeff = {{a0, a1},{0., b1}, poles};
	return coeff;
}

double autoMem_allPass6(double input, double freq)
{
	iir_Coeffs coeffs = Coeffs_allPass6(freq);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		ALL PASS 12

iir_Coeffs Coeffs_allPass12(double freq, double q)
{
	const unsigned int poles = 2;

	double freq1 = freq <= 1. ? freq/2. : freq/iir_sampleRate;
	double theta = PI*freq1;
	double alpha = (tan(theta/q) - 1.)/(tan(theta/q) + 1.);
	double beta  = -cos(2.*theta);
	double a0 = -alpha;
	double a1 = beta*(1. - alpha);
	double a2 = 1.;
	double b1 = -a1;
	double b2 = alpha;

	iir_Coeffs coeff = {{a0, a1, a2},{0., b1, b2}, poles};
	return coeff;
}

double autoMem_allPass12(double input, double freq, double q)
{
	iir_Coeffs coeffs = Coeffs_allPass12(freq, q);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		LOW SHELVING

iir_Coeffs Coeffs_lowShelving(double freq, double gain)
{
	const unsigned int poles = 1;

	double freq1 = freq <= 1. ? freq/2. : freq/iir_sampleRate;
	double theta = PI*freq1;
	double beta = 4./(1. + gain);
	double delta = beta*tan(theta);
	double gamma = (1. - delta)/(1. + delta);
	double a0 = (1. - gamma)/2.;
	double a1 = a0;
	double b1 = gamma;

	iir_Coeffs coeff = {{a0, a1},{0., b1}, poles};
	return coeff;
}

double autoMem_lowShelving(double input, double freq, double gain_dB)
{
	double gain = pow(10., gain_dB/20.);
	iir_Coeffs coeffs = Coeffs_lowShelving(freq, gain);
	return input + (gain - 1.)*autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		HIGH SHELVING

iir_Coeffs Coeffs_highShelving(double freq, double gain)
{
	const unsigned int poles = 1;

	double freq1 = freq <= 1. ? freq/2. : freq/iir_sampleRate;
	double theta = PI*freq1;
	double beta = (1. + gain)/4.;
	double delta = beta*tan(theta);
	double gamma = (1. - delta)/(1. + delta);
	double a0 = (1. + gamma)/2.;
	double a1 = -a0;
	double b1 = gamma;

	iir_Coeffs coeff = {{a0, a1},{0., b1}, poles};
	return coeff;
}

double autoMem_highShelving(double input, double freq, double gain_dB)
{
	double gain = pow(10., gain_dB/20.);
	iir_Coeffs coeffs = Coeffs_highShelving(freq, gain);
	return input + (gain - 1.)*autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		PEAK

iir_Coeffs Coeffs_peak(double freq, double q, double gain)
{
	const unsigned int poles = 2;

	double freq1 = freq <= 1. ? freq/2. : freq/iir_sampleRate;
	double theta = 2.*PI*freq1;
	double sigma = 4./(1. + gain);
	double beta = .5*(1. - sigma*tan(.5*theta/q))/(1. + sigma*tan(.5*theta/q));
	double gamma = (.5 + beta)*cos(theta);
	double a0 = .5 - beta;
	double a1 = 0.;
	double a2 = -a0;
	double b1 = 2.*gamma;
	double b2 = -2.*beta;

	iir_Coeffs coeff = {{a0, a1, a2},{0., b1, b2}, poles};
	return coeff;
}

double autoMem_peak(double input, double freq, double q, double gain_dB)
{
	double gain = pow(10., gain_dB/20.);
	iir_Coeffs coeffs = Coeffs_peak(freq, q, gain);
	return input + (gain - 1.)*autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		PEAK CONSTANT Q

iir_Coeffs Coeffs_peakConstQ(double freq, double q, double gain)
{
	const unsigned int poles = 2;

	double freq1 = freq <= 1. ? freq/2. : freq/iir_sampleRate;
	double theta = PI*freq1;
	double k = tan(theta);
	double d0 = 1. + k/q + k*k;
	double e0 = 1. + k/(gain*q) + k*k;
	double alpha = 1. + gain*k/q + k*k;
	double beta = 2.*(k*k - 1.);
	double gamma = 1. - gain*k/q + k*k;
	double delta = 1. - k/q + k*k;
	double eta = 1. - k/(gain*q) + k*k;

	int boost = gain >= 1.;
	double a0 = boost ? alpha/d0  : d0/e0;
	double a1 = boost ? beta/d0   : beta/e0;
	double a2 = boost ? gamma/d0  : delta/e0;
	double b1 = boost ? -beta/d0  : -beta/e0;
	double b2 = boost ? -delta/d0 : -eta/e0;

	iir_Coeffs coeff = {{a0, a1, a2},{0., b1, b2}, poles};
	return coeff;
}

double autoMem_peakConstQ(double input, double freq, double q, double gain_dB)
{
	double gain = pow(10., gain_dB/20.);
	iir_Coeffs coeffs = Coeffs_peakConstQ(freq, q, gain);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		LINKWITZ-RILEY LP12

iir_Coeffs Coeffs_linkwitzRileyLp12(double freq)
{
	const unsigned int poles = 2;

	double omega = PI*freq;
	double theta = omega/iir_sampleRate;
	double kappa = omega/tan(theta);
	double delta = kappa*kappa + omega*omega + 2.*kappa*omega;
	double a0 = omega*omega/delta;
	double a1 = 2.*a0;
	double a2 = a0;
	double b1 = (2.*kappa*kappa - 2.*omega*omega)/delta;
	double b2 = (2.*kappa*omega - kappa*kappa - omega*omega)/delta;
	
	iir_Coeffs coeff = {{a0, a1, a2},{0., b1, b2}, poles};
	return coeff;
}

double autoMem_linkwitzRileyLp12(double input, double freq)
{
	iir_Coeffs coeffs = Coeffs_linkwitzRileyLp12(freq);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		LINKWITZ-RILEY HP12

iir_Coeffs Coeffs_linkwitzRileyHp12(double freq)
{
	const unsigned int poles = 2;

	double omega = PI*freq;
	double theta = omega/iir_sampleRate;
	double kappa = omega/tan(theta);
	double delta = kappa*kappa + omega*omega + 2.*kappa*omega;
	double a0 = kappa*kappa/delta;
	double a1 = -2.*a0;
	double a2 = a0;
	double b1 = (2.*kappa*kappa - 2.*omega*omega)/delta;
	double b2 = (2.*kappa*omega - kappa*kappa - omega*omega)/delta;
	
	iir_Coeffs coeff = {{a0, a1, a2},{0., b1, b2}, poles};
	return coeff;
}

double autoMem_linkwitzRileyHp12(double input, double freq)
{
	iir_Coeffs coeffs = Coeffs_linkwitzRileyHp12(freq);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		FAST LP6

iir_Coeffs Coeffs_fastLp6(double normalizedFreq)
{
	const unsigned int poles = 1;

	double a0 = normalizedFreq;
	double b1 = 1. - normalizedFreq;

	iir_Coeffs coeff = {{a0},{0., b1}, poles};
	return coeff;
}

double autoMem_fastLp6(double input, double normalizedFreq)
{
	iir_Coeffs coeffs = Coeffs_fastLp6(normalizedFreq);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		FAST HP6

iir_Coeffs Coeffs_fastHp6(double normalizedFreq)
{
	const unsigned int poles = 1;

	double gamma = 1. - 2.*normalizedFreq;
	double a0 = (1. + gamma)/2.;
	double a1 = -a0;
	double b1 = gamma;

	iir_Coeffs coeff = {{a0, a1},{0., b1}, poles};
	return coeff;
}

double autoMem_fastHp6(double input, double normalizedFreq)
{
	iir_Coeffs coeffs = Coeffs_fastHp6(normalizedFreq);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		FAST LP12

iir_Coeffs Coeffs_fastLp12(double normalizedFreq, double damping)
{
	const unsigned int poles = 2;

	double f1 = 2.*normalizedFreq - 1.;
	double sinfApprox = 1. - f1*f1;
	double dsint = damping*sinfApprox/2.;
	double beta  = .5*(1. - dsint)/(1. + dsint);
	double cosfApprox = 1. - 2.*normalizedFreq;
	double gamma = (.5 + beta)*cosfApprox;
	double delta = .5 + beta - gamma;
	double a0 = delta/2.;
	double a1 = delta;
	double a2 = a0;
	double b1 = 2.*gamma;
	double b2 = -2.*beta;
	
	iir_Coeffs coeff = {{a0, a1, a2},{0., b1, b2}, poles};
	return coeff;
}

double autoMem_fastLp12(double input, double normalizedFreq, double damping)
{
	iir_Coeffs coeffs = Coeffs_fastLp12(normalizedFreq, damping);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------
//		FAST HP12

iir_Coeffs Coeffs_fastHp12(double normalizedFreq, double damping)
{
	const unsigned int poles = 2;
	
	double f1 = 2.*normalizedFreq - 1.;
	double sinfApprox = 1. - f1*f1;
	double dsint = damping*sinfApprox/2.;
	double beta  = .5*(1. - dsint)/(1. + dsint);
	double cosfApprox = 1. - 2.*normalizedFreq;
	double gamma = (.5 + beta)*cosfApprox;
	double sigma = .5 + beta + gamma;
	double a0 = sigma/2.;
	double a1 = -sigma;
	double a2 = a0;
	double b1 = 2.*gamma;
	double b2 = -2.*beta;
	
	iir_Coeffs coeff = {{a0, a1, a2},{0., b1, b2}, poles};
	return coeff;
}

double autoMem_fastHp12(double input, double normalizedFreq, double damping)
{
	iir_Coeffs coeffs = Coeffs_fastHp12(normalizedFreq, damping);
	return autoMem_filter(input, &coeffs, IIR_LINEAR);
}

//-----------------------------------------------------------------------------

#endif // AUTOIIR_IMPLEMENTATION

#ifdef __cplusplus
} // extern "C"
#endif

#endif // AUTOIIR_H

