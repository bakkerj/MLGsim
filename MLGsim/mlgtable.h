#pragma once
#ifndef __INC_MLG_MGLTABLE_HH__
#define __INC_MLG_MGLTABLE_HH__

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <cmath>

#include "locus.h"
#include "mlgtuple.h"
#include "datareader.h"

namespace mlg {

  class MLGTable {
	public:
    // for easier access to the tuple (AlleleType = int, SampleType = string)
		typedef MLGTuple<AlleleType, SampleType> TupleType;
    // ploidy enumeration: As long as enum value == number of alleles per locus, all algorithms should work with it
    enum Ploidy {
      Unspecified = 0,		// unspecified ploidy
      Haploid = 1,				// indicates haploid loci
      Diploid = 2,				// indicates diploid loci.
    };

    // constructor: creates a new MLG table; ploidy == Ploidy (currently only Haploid or Diploid), locusCount== number of loci, simCount == number of simulations
    MLGTable(Ploidy ploidy, size_t locusCount, bool freq, bool pgen, size_t simCount)
			: d_ploidy(ploidy), d_locusCount(locusCount), d_loci(locusCount), d_verbose(1), b_freq(freq), b_pgen(pgen), d_simCount(simCount)
    {}

		~MLGTable()
		{}

    // insert a new MLG tuple into the table
		// NB: each sample should have a unique sample id, otherwise the pSex calculations get screwed
    void insert(TupleType &tuple);

    // returns the ploidy
    Ploidy ploidy() const 
		{ return d_ploidy; }

    // returns the number of loci
    size_t locusCount() const 
		{ return d_locusCount; }

    // returns a set of all unique mlgs
    std::vector<TupleType> const &mlgs() const 
		{ return d_mlgSet; }

    // returns a set of all unique mlgs
    std::vector<TupleType> &mlgs() 
		{ return d_mlgSet; }

    // returns a set of all sample ids (one for each inserted sample with a unique sample id)
    std::vector<SampleType> const &samples() const 
		{ return d_samples; }

    // sets the name of a locus
    void setLocusName(size_t index, std::string const &name)
    {
      assert(index < d_loci.size());
      d_loci[index].setName(name);
    }

    // sets the frequencies of this table to the ones in loci
    void setFrequencies(std::vector<Locus> const &loci)
    {
      assert(loci.size() == d_loci.size());
      for (size_t i = 0; i < d_loci.size(); ++i) {
//				assert(loci[i].name() == d_loci[i].name());
				d_loci[i].clearValues();
				for (Locus::container::const_iterator j = loci[i].alleles().begin(); j != loci[i].alleles().end(); ++j)
					d_loci[i].setAlleleValue((*j).first,(*j).second);
				d_loci[i].setFrequencies(true);
      }
    }

    // clears the internal list of pSex values
    void clearPSexValues()
    { d_pSexValues.clear(); }

    // appends a pSex value to the internal pSex value list
    void addPSexValue(double value)
    {
      std::vector<double>::iterator i = std::lower_bound(d_pSexValues.begin(), d_pSexValues.end(), value);
      d_pSexValues.insert(i, value);
    }

    // returns the calculated clonal richness (R = \frac{G-1}{N-1})
    double clonalRichness1()
    { return (double)(d_mlgSet.size() - 1) / (double)(d_samples.size() - 1); }

    // returns the calculated clonal richness (P_d = \frac{G}{N})
    double clonalRichness2()
    { return (double)d_mlgSet.size() / (double)d_samples.size(); }

    // returns the calculated clonal richness (P_de = \frac{Ge}{N})
    double clonalRichness3()
    { return effectiveMLGCount() / (double)d_samples.size(); }

    // returns the critical value at the specified significance level
    double criticalValue(double significanceLevel) const
    {
      assert(significanceLevel >= 0.0 && significanceLevel <= 1.0);
      // if no simulated pSex are generated everything is significant!
      if (d_pSexValues.empty()) return 1;
      return d_pSexValues[(d_pSexValues.size() - 1) * significanceLevel];
    }

    // returns the significance level of this pSex value
    double calcSignificance(double pSex) const
    {
      std::vector<double>::const_iterator i = std::lower_bound(d_pSexValues.begin(), d_pSexValues.end(), pSex);
      if (i == d_pSexValues.end()) 
				return 0.0;
      double i0 = (i - d_pSexValues.begin());
      return i0 / (double)(d_pSexValues.size());
    }

    // returns the number of unique MLGs (i.e. based on a single sample)
    size_t uniqueMLGCount() const
    {
      size_t count = 0;
      for (std::vector<TupleType>::const_iterator i = d_mlgSet.begin(); i != d_mlgSet.end(); ++i)
				count += ((*i).samples().size() == 1 ? 1 : 0);
      return count;
    }

    // returns the number of non-unique MLGs (i.e. based on multiple samples)
    size_t nonUniqueMLGCount() const
    {
      size_t count = 0;
      for (std::vector<TupleType>::const_iterator i = d_mlgSet.begin(); i != d_mlgSet.end(); ++i)
				count += ((*i).samples().size() == 1 ? 0 : 1);
      return count;
    }

    // returns the effective number of MLGs (G_e = 1 / sum of MLG frequencies^2)
    double effectiveMLGCount() const
    {
      double sum = 0.0;
			double total = (double)d_samples.size();
      for (std::vector<TupleType>::const_iterator i = d_mlgSet.begin(); i != d_mlgSet.end(); ++i)
				sum += pow((((double)((*i).samples().size())) / total), 2);
      return 1.0 / sum;
    }

    size_t critCount(double criticalValue) const
    {
      std::pair<std::vector<double>::const_iterator,	std::vector<double>::const_iterator> r = std::equal_range(d_pSexValues.begin(), d_pSexValues.end(), criticalValue);
      return r.second - r.first;
    }

    double critLevel(double significanceLevel) const;

    double critPercentage(double significanceLevel)
    {
      std::pair<std::vector<double>::const_iterator, std::vector<double>::const_iterator> r = std::equal_range(d_pSexValues.begin(), d_pSexValues.end(), criticalValue(significanceLevel));
      return (double)(r.second - r.first) / d_pSexValues.size() * 100.0;
    }

    // returns the sorted list of psex values
    std::vector<double> const &pSexValues() const 
		{ return d_pSexValues; }

    // returns the number of samples that are heterozygote in the specified locus
    size_t countHeterozygotes(size_t lociIndex) const;

    // returns the observed heterozygosity (Ho)
    double observedHeterozygosity(size_t locusIndex) const
    {
			double Ho;
			b_freq ? Ho = (double)countHeterozygotes(locusIndex) / (double)d_samples.size() : Ho = (double)countHeterozygotes(locusIndex) / (double)d_mlgSet.size();
			return Ho;
    }

    // returns the Fis for the specified locus
    double calcFis(size_t locusIndex) const
    {
      assert(locusIndex < d_loci.size());
      double He = d_loci[locusIndex].calcHe();
//			assert(He > 0.0);
			double Ho = observedHeterozygosity(locusIndex);
			double Fis;
			He > 0.0 ? Fis = (He - observedHeterozygosity(locusIndex)) / He : Fis = 1.0;
			return Fis;
    }

		// checks sum of frequencies == 1.0
    void checkFrequencies();

    // calculates the allele frequencies (should be done only once!)
    void calculateFrequencies();

    // calculates pGen value
    double calcPGen(TupleType const &mlg);

    // calculates all pSex values
    void calculatePSexValues(size_t sampleSize);

    // sets the title name of the table
    void setTitle(std::string const &title)
    { d_title = title; }

    // returns the title name of the table
    std::string const &title() const 
		{ return d_title; }

    // print to stream
    void printFrequencies(std::ostream &stream, char sep = '\t');
    void printMLGTable(std::ostream &stream, char sep = '\t');
    void print(std::ostream &stream, char sep = '\t');
    void printOutputTable(std::ostream &stream, char sep = '\t');

    // returns the pSex value distribution
    std::vector<double> const &pSexValues() 
		{ return d_pSexValues; }

    // returns the loci
    std::vector<Locus> const &loci() const 
		{ return d_loci; }

    // sets the verbose level of errors
    void setVerbose(size_t value) 
		{ d_verbose = value; }

//    static MLGTable fromTableSimulation(DataReader &reader, size_t simCount); // called once from Main()
    static MLGTable fromTableSimulation(DataReader &reader); // called once from Main(); number of sims is a parameter!
    static MLGTable loadFromFile(DataReader &reader); // called once from function fromTableSimulation()
    static MLGTable fromSimulation(size_t sampleSize, Ploidy ploidy, std::vector<Locus> const &loci, bool freq, bool pgen, size_t simCount); // called simCount times from function fromTableSimulation()

  private:
    Ploidy d_ploidy;
    size_t d_locusCount;
    std::vector<Locus> d_loci;
    std::vector<TupleType> d_mlgSet;
    std::vector<SampleType> d_samples;
    std::vector<double> d_pSexValues;
    std::string d_title;
    size_t d_verbose;
		size_t d_simCount;
		bool b_freq;
		bool b_pgen;
  };
};

#endif // __INC_MLG_MGLTABLE_HH__
