#pragma once
#ifndef __INC_MLG_CONFIG_HH__
#define __INC_MLG_CONFIG_HH__

namespace mlg {
  // Allel identifier data type
  typedef int AlleleType;
  // Sample identifier data type (should have less-than and equals operator)
  typedef std::string SampleType;
	// Frequency calculation: SAMPLE == use sample size; MLG == use number of MLGs
	#define SAMPLE true
	#define MLG false
	// Frequency calculation: RAW == use uncorrected data; RR == use corrected (round-robin) data 
	#define RAW true
	#define RR false
	// Pgen calculation: HWE == assuming H-W equilibrium; FIS == corrected for departure from HWE using Fis 
	#define FIS true
	#define HWE false
};

#endif // __INC_MLG_CONFIG_HH__
