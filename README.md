# AutoIIR
Easy to use stb-style IIR filter library for audio processing loops. The goal was to get as close as possible to allow the user just call a filter function in processing loop and focus on the creative side of audio programming. 

## Example program

The following pseudo-code demonstrates processing input with 3 static filters in parallel. 

```c
// Required for ONE implementation file
#define AUTOIIR_IMPLEMENTATION

// Optional compile time options
#define AUTOIIR_DEBUG
#define AUTOIIR_MAX_AUDIO_INPUTS 2

// Include the library AFTER setting compile time options
#include <autoiir.h>

//	...
//
//	framework code
//
//	...

void initialize()
{
	// Optional callback setup
	// Does nothing if 'AUTOIIR_DEBUG' is not defined
	iir_setErrorMessageCallback(myErrorPrintingCallback);
	
	// Required for non-normalized input values
	iir_setSampleRate(44100);
}

void processBlock(double** audioBuffer)
{
	// Processing loop
	for (int c = 0; c < audioInputsCount; c++)
		for (int s = 0; s < samplesToProcess; s++)
		{
			// Required for 'autoMem' filters
			// Reset filter counter to map filter state
			autoMem_initProcessing(c);
			
			double x = audioBuffer[c][s];
			
			// Non-normalized frequency
			double lpOut = lowPass12(x, 2000/*Hz*/, .707);
			
			// Normalized frequency = 0.5 * Nyquist = 10025 Hz
			double bsOut = bandStop12(x, 0.5, 1);
			
			// Non-linear filter
			// All filter stages go trough 'atan()'
			iir_FuncTable fn = initFuncTableWith(atan);
			iir_Coeffs coeff = Coeffs_buttHp12(150/*Hz*/);
			double hpOut = autoMem_filter(x, &coeff, &fn);
			
			audioBuffer[c][s] = lpOut + bsOut + hpOut;
		}
}

void reset()
{
	// Reset filter state
	autoMem_clear();
}
```
