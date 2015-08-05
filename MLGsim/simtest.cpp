#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
// added to check OS: <unistd.h> is unix-specific
#ifdef _MSC_VER
  #include <windows.h>
#else
  #include <unistd.h>
#endif
#include "mlgtable.h"

using namespace std;
using namespace mlg;

bool exists(string const &path)
{
  struct stat statBuffer;
  if (!stat(path.c_str(), &statBuffer))
    return true;
  return false;
}

int main(int argc, char *argv[])
{
  for (char **i = argv + 1; *i; ++i) {
		// read & sort input
    ifstream stream(*i);
		if (stream.fail())	// exit met code 2 als file niet (goed) geopend is
			exit(2);
		cout << "Reading file " << *i << endl; //	ok: *i == rrdata.txt
		DataReader reader(stream);
		// simulate 1000 times
//		MLGTable table = MLGTable::fromTableSimulation(reader, 1000);
		MLGTable table = MLGTable::fromTableSimulation(reader);	// number of sims == param in input file
		// write output to file
    string outputFilename = "simresults_";
    string fullTableFilename = "fulltable_";
    if (table.title().empty()) {
      outputFilename += *i;
      fullTableFilename += *i;
    } 
		else {
      outputFilename += table.title();
      fullTableFilename += table.title();
    }
    outputFilename += ".csv";
    fullTableFilename += ".csv";
    if (exists(outputFilename))
      cout << "Warning: overwriting output file '" << outputFilename << "'." << endl;
    if (exists(fullTableFilename))
      cout << "Warning: overwriting output file '" << fullTableFilename << "'." << endl;
		// write results
    ofstream res(outputFilename.c_str());
    res << "Datafile," << *i << endl << endl;
    table.print(res,',');
    res << endl << endl;
    res << "Specified parameters:" << endl;
    for (map<string, string>::const_iterator j = reader.parameters().begin(); j != reader.parameters().end(); ++j)
      res << (*j).first << "," << (*j).second << endl;
    res << endl;
		// write input data
    ofstream full(fullTableFilename.c_str());
    full << "Datafile," << *i << endl << endl;
    table.printOutputTable(full, ',');
    full << endl << endl;
		cout << "Results are written to '" << outputFilename << "'." << endl;
  }
  return 0;
}
