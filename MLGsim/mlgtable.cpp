#include "mlgtable.h"
#include "math.h"
#include <cmath>
#include <queue>
#include <sstream>
#include <iomanip>

using namespace std;
using namespace mlg;

typedef signed int ssize_t;

// Insert tuple into the table
void MLGTable::insert(TupleType &tuple)
{
  // Search for an MLG with the same allele values (binary search in sorted data, so this is fast)
  std::vector<TupleType>::iterator i = std::lower_bound(d_mlgSet.begin(),d_mlgSet.end(),tuple);
  // If there is no MLG with the same alleles, insert a new one, otherwise add the sample id(s) to the existing MLG
	bool added = false;
  if (i != d_mlgSet.end() && (*i) == tuple)
    (*i).insertSampleFrom(tuple);
  else {
    d_mlgSet.insert(i, tuple);
		added = true;
  }
  // Keep the count of all alleles inside a locus, so we can calculate the frequencies
	if (b_freq) { // use default option SAMPLE to calculate allele frequencies from sample size
		for (size_t j = 0; j < d_locusCount; ++j) {
			for (ssize_t k = 0; k < d_ploidy; ++k) {
				d_loci[j].insertAllele(tuple.get(j, k));
			}
		}
	}
	else { // use alternate option MLG to calculate allele frequencies from number of MLGs
		if (added){
			for (size_t j = 0; j < d_locusCount; ++j) {
				for (ssize_t k = 0; k < d_ploidy; ++k) {
					d_loci[j].insertAllele(tuple.get(j, k));
				}
			}
		}
	}
  // Remember all sample ids, so we know how many samples there were
  for (vector<SampleType>::const_iterator s = tuple.samples().begin(); s != tuple.samples().end(); ++s) {
    d_samples.push_back(*s);
  }
}

void MLGTable::printMLGTable(std::ostream &stream, char sep)
{
  // Print header
  stream << "MLG" << sep << "n" << sep << "PSex" << sep << "Significance" << sep << "Level" << sep << "PGen";
  for (std::vector<Locus>::const_iterator i = d_loci.begin(); i != d_loci.end(); ++i)
    for (ssize_t j = 0; j < d_ploidy; ++j)
      stream << sep << (*i).name() << "_" << j;
  stream << std::endl;
  size_t mlgId = 1;
  for (std::vector<TupleType>::const_iterator i = d_mlgSet.begin(); i != d_mlgSet.end(); ++i) {
    stream << "MLG" << mlgId++ << sep << (*i).samples().size();
    if ((*i).samples().size() > 1) {
      string slev;
      if ((*i).pSex() < criticalValue(0.05)) {
				slev += "*";
				if ((*i).pSex() < criticalValue(0.01)) {
					slev += "*";
					if ((*i).pSex() < criticalValue(0.001)) {
						slev += "*";
					}
				}
      } 
			else 
				slev = "ns";
      stream << sep << (*i).pSex() << sep << calcSignificance((*i).pSex()) << sep << slev	<< sep << (*i).pGen();
    } 
		else
      stream << sep << "\"-\"" << sep << "\"-\"" << sep << "\"-\"" << sep << (*i).pGen();
    for (TupleType::const_iterator j = (*i).begin(); j != (*i).end(); ++j)
      stream << sep << *j;
    stream << std::endl;
  }
  // Print locus info
  stream << endl << sep << sep << sep;
  for (std::vector<Locus>::const_iterator i = d_loci.begin(); i != d_loci.end(); ++i)
    stream << sep << (*i).name() << sep;
  stream << std::endl;
  stream << "Number of alleles" << sep << sep << sep;
  for (size_t i = 0; i < d_loci.size(); ++i)
    stream << sep << d_loci[i].uniqueAlleleCount() << sep;
  stream << std::endl;
  stream << "He" << sep << sep << sep;
  for (size_t i = 0; i < d_loci.size(); ++i)
    stream << sep << d_loci[i].calcHe() << sep;
  stream << std::endl;
  stream << "Ho" << sep << sep << sep;
  for (size_t i = 0; i < d_loci.size(); ++i)
	  stream << sep << observedHeterozygosity(i) << sep;
  stream << std::endl;
  stream << "Fis" << sep << sep << sep;
  for (size_t i = 0; i < d_loci.size(); ++i)
	  stream << sep << calcFis(i) << sep;
  stream << std::endl;
}

void MLGTable::printOutputTable(std::ostream &stream, char sep)
{
  // Print header
  stream << "Sample" << sep << "MLG" << sep << "PSex" << sep << "PValue";
  for (std::vector<Locus>::const_iterator i = d_loci.begin(); i != d_loci.end(); ++i)
    for (ssize_t j = 0; j < d_ploidy; ++j)
      stream << sep << (*i).name() << "_" << j;
  stream << std::endl;
  size_t mlgId = 1;
  for (std::vector<TupleType>::const_iterator i = d_mlgSet.begin(); i != d_mlgSet.end(); ++i) {
    for (std::vector<SampleType>::const_iterator j = (*i).samples().begin(); j != (*i).samples().end(); ++j) {
      stream << *j << sep << "MLG_" << mlgId;
      if ((*i).samples().size() > 1)
				stream << sep << (*i).pSex();
      else
				stream << sep << "\"-\"";
      stream << sep << (*i).pValue();
      for (TupleType::const_iterator k = (*i).begin(); k != (*i).end(); ++k)
				stream << sep << *k;
      stream << std::endl;
    }
    ++mlgId;
  }
}


void MLGTable::printFrequencies(std::ostream &stream, char sep)
{
  stream << endl << "Allele frequencies:" << endl;
  for (std::vector<Locus>::const_iterator i = d_loci.begin(); i != d_loci.end(); ++i) {
		stream << "Locus " << (*i).name() << endl;
    for (std::map<AlleleType, double>::const_iterator j = (*i).alleles().begin(); j != (*i).alleles().end(); ++j)
      stream << sep << (*j).first << sep << (*j).second << endl;
  }
}

void MLGTable::print(std::ostream &stream, char sep)
{
  printMLGTable(stream, sep);
  stream << endl << endl
	 << "Clonal Richness" << endl
	 << "G" << sep << d_mlgSet.size() << endl
	 << "N" << sep << d_samples.size() << endl
	 << "R" << sep << clonalRichness1() << endl
	 << "Pd" << sep << clonalRichness2() << endl
	 << "Ge" << sep << effectiveMLGCount() << endl
	 << "Pde" << sep << clonalRichness3() << endl
	 << "Single MLGs" << sep << uniqueMLGCount() << endl
	 << "Clonal MLGs" << sep << nonUniqueMLGCount() << endl
	 << endl;
  stream << "PValue(0.01)" << sep << criticalValue(0.01) << sep  << "significance level" << sep << critLevel(0.01) << endl;
  stream << "PValue(0.05)" << sep << criticalValue(0.05) << sep  << "significance level" << sep << critLevel(0.05) << endl;
  stream << "PValue(0.10)" << sep << criticalValue(0.10) << sep  << "significance level" << sep << critLevel(0.10) << endl << endl;
  stream << "Samplesize" << sep << d_samples.size() << endl;
  stream << "PValues" << sep << d_pSexValues.size() << endl << endl;
  printFrequencies(stream, sep);
  stream << endl << endl << "PSex Values:" << endl;
  for (std::vector<double>::const_iterator i = d_pSexValues.begin(); i != d_pSexValues.end(); ++i)
    stream << *i << sep << endl;
  stream << endl << endl;
}

void MLGTable::checkFrequencies()
{
  bool error = false;
  for (vector<Locus>::iterator i = d_loci.begin(); i != d_loci.end(); ++i) {
    double sum = 0.0;
    for (Locus::container::const_iterator a = (*i).alleles().begin(); a != (*i).alleles().end(); ++a)
      sum += (*a).second;
    if (fabs(sum - 1.0) > 0.00000001) {
      error = true;
      if (d_verbose > 2)
				cerr << "error: frequency sums up to " << sum << "." << endl;
    }
  }
  if (error)
    cerr << "error: not all frequencies sum up to +-1.0." << endl;
}

void MLGTable::calculateFrequencies()
{
  // Make frequencies from counts
  for (vector<Locus>::iterator i = d_loci.begin(); i != d_loci.end(); ++i) {
    if ((*i).frequencies())
			break;
		b_freq ? (*i).calculateFrequencies(d_samples.size(), d_ploidy) : (*i).calculateFrequencies(d_mlgSet.size(), d_ploidy);
  }
  checkFrequencies();
}

size_t MLGTable::countHeterozygotes(size_t locusIndex) const
{
  size_t count = 0;
  for (std::vector<TupleType>::const_iterator s = d_mlgSet.begin(); s != d_mlgSet.end(); ++s)
		count += (*s).isHeterozygote(locusIndex) ? (b_freq ? (*s).samples().size() : 1) : 0;
  return count;
}

// Simple Pgen calculation (no round-robin freqs)
// NB: b_pgen = TRUE: use Fis to adjust for departures from H-W equilibrium, otherwise assume H-W equilibrium
double MLGTable::calcPGen(TupleType const &mlg)
{
  double p = 1.0, freq;
	for (size_t i = 0; i < d_loci.size(); ++i) {
		freq = 1.0;
		for (ssize_t a = 0; a < d_ploidy; ++a)
	     freq *= d_loci[i].allele(mlg.get(i, a));
		if (b_pgen) {
			double fis = calcFis(i);
			(mlg.isHeterozygote(i)) ? freq *= (1.0 - fis) : freq = freq * (1.0 - fis) + d_loci[i].allele(mlg.get(i, 0)) * fis;
		}
		if (mlg.isHeterozygote(i))
			freq *= 2.0;
		p *= freq;
  }
  return p;
}

void MLGTable::calculatePSexValues(size_t sampleSize)
{
  calculateFrequencies();
  for (std::vector<TupleType>::iterator s = d_mlgSet.begin(); s != d_mlgSet.end(); ++s) {
    if ((*s).samples().size() < 1) { // Only calculate for MLG with multiple samples
      (*s).setPSex(0.0);
      continue;
    }
    double res = 0.0;
    double pGen = calcPGen(*s);
    // For all samples
    for (ssize_t i = (*s).samples().size() - 1; i >= 0; --i) {
      res += Math::binom(sampleSize, i) * pow(pGen, double(i)) * pow(1.0 - pGen, double(sampleSize - i));
    }
    // Due to lack of precision in the input data or, for a large data set, lack of precision in the calculations,
    // tiny errors can occur in the resulting probability that might cause it to get larger than 1.0 
    if (res > 1.0) {
			if (res - 1.0 > 0.0000001) {
				cerr << "Warning: rounded " << setprecision(12) << res << " to 1.0." << endl;
//				cout << "res: " << setprecision(12) << res << endl;
			}
      res = 1.0;
    }
    (*s).setPGen(pGen);
    (*s).setPSex(1.0 - res);
    d_pSexValues.push_back(1.0 - res);
  }
  sort(d_pSexValues.begin(), d_pSexValues.end());
}

MLGTable MLGTable::loadFromFile(DataReader &reader)
{
  std::string locusNames;
  std::string title;
  Ploidy ploidy = Unspecified;
//  Ploidy ploidy = Diploid;  // gebruikt om zonder haplo optie te runnen; program werkt dan niet, maar Pgen is niet gedef voor haplo
	bool pgen = HWE, freq = SAMPLE;
  size_t locusCount = 0, simCount = 0;
  if (reader.hasParameter("PLOIDY")) {
    string value = reader.parameter("PLOIDY");
    if (value == "DIPLOID" || value == "2")
      ploidy = Diploid;
    else if (value == "HAPLOID" || value == "1")
      ploidy = Haploid;
    else
      std::cerr << "error: only HAPLOID and DIPLOID are supported for now (not '" << value << "')." << std::endl;
  }
	if (reader.hasParameter("FREQUENCY")) {
		string value = reader.parameter("FREQUENCY");
    if (value == "SAMPLE")
      freq = true;
    else if (value == "MLG")
      freq = false;
    else
      std::cerr << "error: FREQUENCY is unspecified, assuming SAMPLE." << std::endl;
  }
	if (reader.hasParameter("MODEL")) {
		string value = reader.parameter("MODEL");
    if (value == "FIS")
      pgen = true;
    else if (value == "HWE")
      pgen = false;
    else
      std::cerr << "error: MODEL is unspecified, assuming HWE." << std::endl;
  }
  if (reader.hasParameter("SIMULATIONS"))
      simCount = std::atoi(reader.parameter("SIMULATIONS").c_str());
  if (reader.hasParameter("LOCUSCOUNT"))
      locusCount = std::atoi(reader.parameter("LOCUSCOUNT").c_str());
  if (reader.hasParameter("LOCUSNAMES"))
    locusNames = reader.parameter("LOCUSNAMES") + ",";
  if (reader.hasParameter("TITLE"))
    title = reader.parameter("TITLE");
  // Check input
  if (ploidy == Unspecified) {
    std::cerr << "error: PLOIDY is unspecified, assuming DIPLOID." << std::endl;
    ploidy = Diploid;
  } 
  if (locusCount == 0) {
    std::cerr << "error: LOCUSCOUNT is unspecified, assuming 1." << std::endl;
    locusCount = 1;
  }
  if (simCount == 0) {
    std::cerr << "error: SIMULATIONS is unspecified, assuming 1000." << std::endl;
    simCount = 1000;
  }
  // Initialize MLGTable
  MLGTable res(ploidy, locusCount, freq, pgen, simCount);
  res.setTitle(title);
  // Parse locus names
  size_t idx = 0;
  for (size_t i = 0, e = locusNames.find(','); e != std::string::npos; i = e + 1, e = locusNames.find(',', e + 1), ++idx) {
    if (idx >= res.d_loci.size())
      std::cerr << "error: too many locus names given in LOCUSNAMES." << std::endl;
    else
      res.setLocusName(idx, locusNames.substr(i, e - i));
  }
  SampleType sample;
  if (reader.columnCount() != ploidy*locusCount + 1) {
    cerr << "error: not enough columns in table according to ploidy and locus count (should be: " << (ploidy * locusCount + 1) << " has: " << reader.columnCount() << ")." << endl;
    return res;
  }
  // Read data
  while (reader.next()) {
    MLGTuple<AlleleType, SampleType> tuple(ploidy, locusCount);
    tuple.insertSample(reader.get(0));
    for (size_t alleleI = 0; alleleI < ploidy * locusCount; ++alleleI)
      tuple[alleleI] = atoi(reader.get(alleleI + 1).c_str());
    res.insert(tuple);
  }
  res.calculatePSexValues(res.samples().size());
  return res;
}

MLGTable MLGTable::fromSimulation(size_t sampleSize, Ploidy ploidy, vector<Locus> const &loci, bool freq, bool pgen, size_t simCount)
{
  MLGTable res(ploidy, loci.size(), freq, pgen, simCount);
  for (vector<Locus>::const_iterator l = loci.begin(); l != loci.end(); ++l)
    res.setLocusName(l - loci.begin(), (*l).name());
  // single simulation: randomly reallocate alleles per locus for each sample to create a new MLG table
  for (size_t s = 0; s < sampleSize; ++s) {
    TupleType tuple(ploidy,loci.size());
    stringstream ss;
    ss << "s" << s;
    tuple.insertSample(ss.str());
    for (vector<Locus>::const_iterator l = loci.begin(); l != loci.end(); ++l) {
      priority_queue<AlleleType> alleles;
      for (ssize_t p = 0; p < ploidy; ++p) {
				double r = Math::fRandom();
				for (std::map<AlleleType,double>::const_iterator a = (*l).alleles().begin(); a != (*l).alleles().end(); ++a) {
					if (r < (*a).second) {
						alleles.push((*a).first);
						break;
					} 
					else
						r -= (*a).second;
				}
      }
      // Put the alleles sorted from low to high in the tuple
      for (;!alleles.empty(); alleles.pop())
				tuple.set(l - loci.begin(), alleles.size() - 1, alleles.top());
    }
    res.insert(tuple);
  }
  res.calculatePSexValues(res.samples().size());
  return res;
}

//MLGTable MLGTable::fromTableSimulation(DataReader &reader, size_t simCount)
MLGTable MLGTable::fromTableSimulation(DataReader &reader)
{
  // load input
  MLGTable input = MLGTable::loadFromFile(reader);
  // reset input
  input.clearPSexValues();
	// simulate simCount times:
	size_t simCount = input.d_simCount;
  for (size_t c = 0; c < simCount; ++c) {
//		if (c % 10 == 0)
//			cout << "simulating run " << c << endl;
    MLGTable simulated = fromSimulation(input.samples().size(), input.ploidy(), input.loci(), input.b_freq, input.b_pgen, simCount);
    for (std::vector<TupleType>::iterator i = simulated.mlgs().begin(); i != simulated.mlgs().end(); ++i) {
      if ((*i).samples().size() <= 1) // skip MLG based on single sample?
				continue;
      MLGTable table(input.ploidy(),input.loci().size(), input.b_freq, input.b_pgen, simCount);
      table.insert(*i);
      table.setFrequencies(input.loci());
      table.calculatePSexValues(input.samples().size());
      input.addPSexValue(table.mlgs().front().pSex());
    }
  }
  return input;
}

double MLGTable::critLevel(double significanceLevel) const
{
  if (d_pSexValues.empty())
		// there are no datapoints, so 0% lies below the threshold.
    return 0.0;
	return (Math::round(d_pSexValues.size() * significanceLevel) - 1.0) / (double)d_pSexValues.size();
}
