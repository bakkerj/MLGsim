#pragma once
#ifndef __INC_MGL_DATAREADER_HH__
#define __INC_MGL_DATAREADER_HH__

#include <iosfwd>
#include <map>
#include <vector>
#include <string>
#include <stdexcept>
#include <cassert>

namespace mlg {

  class DataReader {
  public:
    // constructor: initializes a reader for the specified stream
    DataReader(std::istream &stream) 
			: d_stream(stream), d_lineNo(1), d_verbose(1)
    { parseHeader(); }

    // moves on to the next data line
    bool next();

    // returns true if the data has a column header
    bool hasColumnHeader() const 
		{ return !d_columnHeader.empty(); }

    // returns the current line number
    size_t lineNo() const 
		{ return d_lineNo; }

    // returns the data of the current line for the specified column index
    std::string const &get(size_t columnIndex) const
    {
      assert(columnIndex < d_columnData.size());
      return d_columnData[columnIndex];
    }

    // returns the data of the current line for the specified column name
    std::string const &get(std::string const &columnName) const
    {
      std::map<std::string,size_t>::const_iterator i = d_columnHeader.find(columnName);
			if (i == d_columnHeader.end())
				throw std::runtime_error("unknown column name '" + columnName + "'");
      return get((*i).second);
    }

    // returns the column header map
    std::map<std::string,size_t> const &columnHeader() const 
		{ return d_columnHeader; }

    // returns true if the given parameter is specified in the stream
    bool hasParameter(std::string const &name)
    { return d_paramHeader.find(name) != d_paramHeader.end(); }

    // returns the parameter with the specified name
    std::string parameter(std::string const &name)
    {
      assert(hasParameter(name));
      return d_paramHeader[name];
    }

    // returns the parameters map
    std::map<std::string,std::string> const &parameters() const
    { return d_paramHeader; }

    // returns the number of columns in the header
    size_t columnCount() const 
		{ return d_columnHeader.size(); }

    // sets the verbose level of errors
    void setVerbose(size_t value) 
		{ d_verbose = value; }

  private:
    void parseHeader();
    bool parseLine();
    std::istream &d_stream;
    std::map<std::string, std::string> d_paramHeader;
    std::map<std::string, size_t> d_columnHeader;
    std::vector<std::string> d_columnData;
    size_t d_lineNo;
    size_t d_verbose;
  };
};

#endif // __INC_MGL_DATAREADER_HH__
