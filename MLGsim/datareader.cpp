#include "datareader.h"
#include <iostream>
#include <cstdio>

using namespace std;
using namespace mlg;

bool DataReader::next()
{
  return parseLine();
}

void DataReader::parseHeader()
{
  string line;
  // Parse parameter header
  while (true) {
    int c = d_stream.get();
    if (c == '#') {
      getline(d_stream, line); // Skip comment line
      ++d_lineNo;
      continue;
    } 
		else if (c == '\r') {
      continue; // Skip line feed
    } 
		else if (c == '\n') {
      ++d_lineNo;
      break; // empty line
    } 
		else if (c == '\t' || c == ' ') {
			if (d_verbose > 2)
				cout << "Warning:" << d_lineNo << ": found white-space character at the beginning of the line, assuming end of header." << endl;
      getline(d_stream, line); // Skip the rest of the line
      break;
    } 
		else {
      string key, value;
      // Parse key.
      for (;c != '='; c = d_stream.get()) {
				if (c == '\n') {
					++d_lineNo;
					cout << "Error:" << d_lineNo << ": header variable with no '=' (line = '" << key << "')." << endl;
					break;
				}
				key += c;
			}
      // Parse value
      for (c = d_stream.get(); c != '\n'; c = d_stream.get()) {
				if (c == EOF) {
					cout << "Error:" << d_lineNo << ": end of file while parsing header variable." << endl;
					return;
				} 
				else if (c == '\t') {
					if (d_verbose > 2)
						cout << "Warning:" << d_lineNo << ": found a tab character in the header variable value field, skipping remainder of the line. If this is an excel artifact, then there is probably nothing to worry about." << endl;
					getline(d_stream, line);
					break;
				}
				value += c;
      }
      ++d_lineNo;
      // Check if it already exists
      if (d_paramHeader.find(key) != d_paramHeader.end())
				cout << "error:" << d_lineNo << ": duplicate key in header ('" << key << "')." << endl;
      d_paramHeader[key] = value;
    }
  }
  size_t columnCount = 0;
  string columnName;
  // Skip comment lines
  while (d_stream.peek() == '#') 
		getline(d_stream, line);
  // Parse data header
  for (int c = d_stream.get(); c != EOF; c = d_stream.get()) {
		if (c == '\r')
			continue;
    // Data is comma or tab separated
		if (c == ',' || c == '\t' || c == '\n') {
      if (!columnName.empty()) {
				if (d_columnHeader.find(columnName) != d_columnHeader.end())
					cout << "error:" << d_lineNo << ": duplicate column name '" << columnName << "'." << endl;
				d_columnHeader[columnName] = columnCount++;
      }
      // End of line
      if (c == '\n') {
				++d_lineNo;
				break;
      }
      columnName.clear();
    } 
		else
      columnName += c;
  }
  // Set column size
  d_columnData.resize(columnCount);
}

bool DataReader::parseLine()
{
  size_t columnIndex = 0;
  string columnData;
  if (!d_stream || d_stream.peek() == EOF)
    return false;
  if (d_stream.peek() == '\t' || d_stream.peek() == ' ') {
    if (d_verbose > 2)
      cout << "warning:" << d_lineNo << ": line starts with a white-space, assuming end of data." << endl;
    return false;
  }
  // Parse data header
  for (int c = d_stream.get();; c = d_stream.get()) {
    if (c == '\r')
			continue;
    // Data is comma or tab separated
    if (c == ',' || c == '\t' || c == '\n' || c == EOF) {
      if (columnIndex < d_columnData.size()) {
				d_columnData[columnIndex++] = columnData;
				columnData.clear();
      } 
			else if (c == EOF) {
				break;
      } 
			else if (!columnData.empty()) {
				cout << "error:" << d_lineNo << ": data line has too many entries." << endl;
				break;
      }
      if (c == '\n') {
				++d_lineNo;
				break;
      }
    } 
		else
      columnData += c;
  }
  if (columnIndex == 0)
		return false;
  else if (columnIndex < d_columnData.size()) {
    cout << "error:" << d_lineNo << ": data line has not enough entries (" << columnIndex << " of " << d_columnData.size() << ")." << endl;
    return false;
  }
  return true;
}
