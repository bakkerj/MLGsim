#pragma once
#ifndef __INC_MLG_MGLTUPLE_HH__
#define __INC_MLG_MGLTUPLE_HH__

#include <cassert>
#include <vector>
#include <set>

namespace mlg {

	// defined in config.h: AlleleType = int, SampleType = string 
	template <typename AlleleType, typename SampleType>
  class MLGTuple {
  public:
    // constructor: ploidy == number of alleles per locus, locusCount == number of loci
    MLGTuple(size_t ploidy, size_t locusCount)
			: d_ploidy(ploidy), d_alleles(ploidy * locusCount), d_pSex(0.0), d_pValue(0.0)
		{}

    // returns true if the alleles of this mlg are the same as those of other
    bool operator==(MLGTuple const &other) const
    {
      assert(other.d_alleles.size() == d_alleles.size());
      return equal(d_alleles.begin(), d_alleles.end(), other.d_alleles.begin());
    }

    // returns true if the alleles of this mlg sort lower than those of other
    bool operator<(MLGTuple const &other) const
    {
      assert(other.d_alleles.size() == d_alleles.size());
      return lexicographical_compare(d_alleles.begin(), d_alleles.end(), other.d_alleles.begin(), other.d_alleles.end());
    }
    
    // returns the allele with the specified index (which should be locusIndex * ploidy + alleleIndex)
    AlleleType operator[](size_t index) const
    {
      assert(index < d_alleles.size());
      return d_alleles[index];
    }

    // returns the allele with the specified index (which should be locusIndex * ploidy + alleleIndex)
    AlleleType &operator[](size_t index)
    {
      assert(index < d_alleles.size());
      return d_alleles[index];
    }

    // returns allele number alleleIndex from locus locusIndex
    AlleleType get(size_t locusIndex, size_t alleleIndex) const
    {
      assert(alleleIndex < d_ploidy);
      size_t index = locusIndex * d_ploidy + alleleIndex;
      assert(index < d_alleles.size());
      return d_alleles[index];
    }

    // sets allele number alleleIndex from locus locusIndex to value
    void set(size_t locusIndex, size_t alleleIndex, AlleleType value)
    {
      assert(alleleIndex < d_ploidy);
      size_t index = locusIndex * d_ploidy + alleleIndex;
      assert(index < d_alleles.size());
      d_alleles[index] = value;
    }

    // inserts a sample to this mlg (doubles are not added)
    void insertSample(SampleType const &sample)
    { d_samples.push_back(sample); }

    // inserts all samples from another MLG tuple
    void insertSampleFrom(MLGTuple const &tuple)
    { d_samples.insert(d_samples.end(), tuple.samples().begin(), tuple.samples().end()); }

    // returns the samples associated with this MLG
    std::vector<SampleType> const &samples() const 
		{ return d_samples; }
/*
    // returns the genotype class associated with this MLG
    size_t genoTypeClass() const 
		{ return d_genoTypeClass; }

    // sets the genotype class associated with this MLG
    void setGenoTypeClass(size_t value) 
		{ d_genoTypeClass = value; }
*/
    // returns the psex value
    double pSex() const 
		{ return d_pSex; }

    // returns the pgen value
    double pGen() const 
		{ return d_pGen; }

    // sets the psex value
    void setPSex(double value) 
		{ d_pSex = value; }

    // sets the pgen value
    void setPGen(double value) 
		{ d_pGen = value; }

    // sets the P value
    void setPValue(double value) 
		{ d_pValue = value; }

    // returns the P value, if calculated
    double pValue() const 
		{ return d_pValue; }

    // returns true if the specified locus is heterozygous
    bool isHeterozygote(size_t locusIndex) const
    { return (d_ploidy == 2 && (get(locusIndex, 0) != get(locusIndex, 1))); }

    typedef typename std::vector<AlleleType>::const_iterator const_iterator;
    
		// returns an iterator to the beginning of the allele list
    const_iterator begin() const 
		{ return d_alleles.begin(); }
   
		// returns an iterator to the end of the allel list
    const_iterator end() const 
		{ return d_alleles.end(); }

  private:
    size_t d_ploidy; // haploid = 1, diploid = 2, etc.
    std::vector<AlleleType> d_alleles;
    std::vector<SampleType> d_samples;
//    size_t d_genoTypeClass;
    double d_pSex;
    double d_pGen;
    double d_pValue;
  };
};

#endif // __INC_MLG_MGLTUPLE_HH__
