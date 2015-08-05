#pragma once
#ifndef __INC_MLG_LOCUS_HH__
#define __INC_MLG_LOCUS_HH__

#include <cstdlib>
#include <algorithm>
#include <map>
#include <cassert>
#include <iostream>
#include "config.h"

namespace mlg {

	// def in config.h: AlleleType == int; SampleType == std::string 

  class Locus {
  public:
    typedef std::map<AlleleType, double> container;
    
		// constructor
		Locus()
			: d_sampleCount(0), d_freqs(false)
		{}

    // sets the locus name
    void setName(std::string const &name) 
		{ d_name = name; }

    // returns the locus name
    std::string const &name() const 
		{ return d_name; }

    // returns the map of allele counts or frequencies (depending on whether calculateFrequencies has been called)
    container const &alleles() const
    { return d_alleles; }

    // returns allele frequency or count (depending on whether calculateFrequencies has been called) of the specified allele
    double allele(AlleleType allele) const
    {
      container::const_iterator i = d_alleles.find(allele);
      return (i != d_alleles.end() ? (*i).second : 0.0);
    }

    // inserts an allele into the locus (keeps the count of all unique alleles)
    void insertAllele(AlleleType a)
    {
      container::iterator i = d_alleles.find(a);
      if (i != d_alleles.end())
				(*i).second += 1.0;
      else
				d_alleles[a] = 1.0;
      ++d_sampleCount;
    }

    // clears all gathered frequency or count data
    void clearValues()
    { d_alleles.clear(); }

    // sets the value (frequency or count) of the specified allele
    void setAlleleValue(AlleleType a, double value)
    { d_alleles[a] = value; }

    // returns the number of unique alleles.
    size_t uniqueAlleleCount() const 
		{ return d_alleles.size(); }

    // calculates the frequencies of the alleles given the sample size and ploidy
		// NB: P&W do not use actual sample size, but number of MLGs * ploidy to calc frequencies
    void calculateFrequencies(size_t sampleSize, size_t ploidy)
    {
      assert(d_sampleCount == sampleSize * ploidy);
      for (container::iterator j = d_alleles.begin(); j != d_alleles.end(); ++j)
				(*j).second /= (double)(d_sampleCount);
      d_freqs = true;
    }

    // returns expected heterozygosity He
		// Note: calculateFrequencies() should be run beforehand!
    double calcHe() const
    {
      double sum = 0.0;
      for (container::const_iterator i = d_alleles.begin(); i != d_alleles.end(); ++i)
				sum += (*i).second * (*i).second;
      return 1.0 - sum;
    }

		// returns true if calculateFrequencies() has been called
    bool frequencies() const { return d_freqs; }

    // set to true to indicate that the allele values are frequencies
    void setFrequencies(bool value) { d_freqs = value; }

  private:
    std::string d_name;
    container d_alleles;
    size_t d_sampleCount;
    bool d_freqs;
  };
};

#endif // __INC_MLG_MGLTABLE_HH__
