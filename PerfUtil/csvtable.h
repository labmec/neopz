/**
 * @file
 * @brief Contains macros and functions to support a string csv table.
 */
#ifndef CSVTABLE_H

#include<cctype>
#include<iostream>// istream
#include <sstream> //stringstream
#include<vector>  // vector
#include <stdlib.h>


/** Template class to create CSVTable columns where the cell values are of type T. */
template<class T>
class CSVTableColumn {

 public:

  CSVTableColumn(const std::string& header)
  {
    setHeader(header);
  }

  /* Gets the column name. */
  const std::string& getHeader() const
  {
    return header;
  }

  /* Sets the column name. */
  void setHeader(const std::string& h)
  {
    header = h;
  }

  /* Checks whether the header name is equal to h. */
  bool isHeader(const std::string& h) const
  {
    return (header == h);
  }
  
  /* Returns the number of rows. Headers do not count as rows. */
  unsigned nRows() const
  {
    return data.size();
  }
  
  /** Add nrows rows to the column. 
   * Return the total number of rows.
   */
  unsigned addRows(unsigned nrows)
  {
    data.resize(data.size() + nrows);
    return data.size();
  }

  /** Returns a reference to the element. */
  T& operator[](unsigned i)
  {
    if (i >= data.size()) {
      std::cerr << "Error(" __FILE__ ":" << __LINE__ << "): Trying to access row " << i
	   << " but column has only " << data.size() << " rows." << std::endl;
      exit(1);
    }
    return data[i];
  }

  /** Returns a constant reference to the element. */
  const T& operator[](unsigned i) const
  {
    if (i >= data.size()) {
      std::cerr << "Error(" __FILE__ ":" << __LINE__ << "): Trying to access row " << i
	   << " but column has only " << data.size() << " rows." << std::endl;
      exit(1);
    }
    return data[i];
  }

 protected:

  std::string header;
  std::vector<T> data;
    
};

template<class Tvar>
class CSVTable
{
 public:
 CSVTable() : nrows(0)  {}
  
  /** 
   * Adds a new row. Returns the new row index (idx).  Do not add rows
   * if there are no columns.
   */
  unsigned addRows(unsigned n=1)
  {
    typename std::vector<CSVTableColumn<Tvar> >::iterator it;
    nrows += n;
    for (it = columns.begin(); it != columns.end(); it++) {
      CSVTableColumn<Tvar>& col = *it;
      col.addRows(n);
      //assert col.nRows() == nrows
    }
    return (nrows-1);
  }

  /** Returns the index of the column named "column_name". */
  int getColIdx(const std::string& column_name)
  {
    int colidx = 0;
    typename std::vector<CSVTableColumn<Tvar> >::iterator it;
    for (it = columns.begin(); it != columns.end(); it++, colidx++) {
      if (it->isHeader(column_name))
	return colidx;
    }
    return -1;
  }

  unsigned addColumn(const std::string& col_name)
  {
    CSVTableColumn<Tvar> col(col_name);
    columns.push_back(col);
    columns.back().addRows(nrows);
    return (columns.size()-1);
  }
  
  /**
   * Sets the value of the cell at column col_id and row row_id.
   * Return 0 if ok, != 0 if error. 
   * Error code: -1 : Invalid row
   *             -2 : Invalid column and createNewCol == false
   */
  int setCell(unsigned row_idx, std::string col_name, std::string value, bool createNewCol = false)
  {
    if (nrows <= row_idx) return -1; // Return -1 row idx is out of bounds.

    int col_idx = getColIdx(col_name);

    if (col_idx < 0) {
      if(!createNewCol) 
	return col_idx; // Return if error.
      else {
	col_idx = addColumn(col_name); 
      }
    }
    CSVTableColumn<Tvar>& col = columns[col_idx];
    col[row_idx] = value;
    return 0;
  }
  int setCell(unsigned row_idx, std::string col_name, int value, bool createNewCol = false)
  {
    std::stringstream str; str << value;
    return setCell(row_idx, col_name, str.str(), createNewCol);
  }
  int setCell(unsigned row_idx, std::string col_name, long value, bool createNewCol = false)
  {
    std::stringstream str; str << value;
    return setCell(row_idx, col_name, str.str(), createNewCol);
  }
  int setCell(unsigned row_idx, std::string col_name, bool value, bool createNewCol = false)
  {
    std::stringstream str; str << value;
    return setCell(row_idx, col_name, str.str(), createNewCol);
  }
  int setCell(unsigned row_idx, std::string col_name, unsigned value, bool createNewCol = false)
  {
    std::stringstream str; str << value;
    return setCell(row_idx, col_name, str.str(), createNewCol);
  }
  int setCell(unsigned row_idx, std::string col_name, double value, bool createNewCol = false)
  {
    std::stringstream str; str.precision(6); str << std::fixed << value;
    return setCell(row_idx, col_name, str.str(), createNewCol);
  }
  int setCell(unsigned row_idx, std::string col_name, float value, bool createNewCol = false)
  {
    std::stringstream str; str << value;
    return setCell(row_idx, col_name, str.str(), createNewCol);
  }

  unsigned nRows() const
  {
    return nrows;
  }

  void write(std::ostream& os) const
  {
    const char* comma = "";
    typename std::vector<CSVTableColumn<Tvar> >::const_iterator it;

    for (it = columns.begin(); it != columns.end(); it++) {
      os << comma << it->getHeader();
      comma = ",";
    }

    for (unsigned i=0; i < nrows; i++) {
      comma = "";
      os << std::endl;
      for (it = columns.begin(); it != columns.end(); it++) {
	os << comma << (*it)[i];
	comma = ",";
      }
    }
  }  

  /** Reads the table from the input stream. 
   *  Returns: 0  if ok
   *           -1 if could not read headers.
   *           1  if number of cells in one or more rows does not match the number of column headers
   */
  int read(std::istream& is)
  {
    std::string line;
    std::vector<std::string> tokens;
    int ret = 0;

    /* reset table */
    columns.clear();
    nrows = 0;

    if (is.eof()) return -1;

    /* read the column headers. */
    getline(is, line);
    if (line.size() == 0) return -1;
    
    split(tokens, line, ",");
    if(tokens.size() == 0) return -1;
    
    /* Build columns. */ 
    std::vector<std::string>::iterator tit;
    for (tit = tokens.begin(); tit != tokens.end(); tit++) {
      std::string hd = *tit;
      trim(hd);
      addColumn(hd);
    }

    unsigned ncols = columns.size();

    while(!is.eof()) {
      std::string line;
      std::vector<std::string> tokens;
      /* read the row. */
      getline(is, line);
      //cout << "Reading line " << nrows << ": " << line << endl;
      split(tokens, line, ",");
      
      unsigned ntokens = tokens.size();
      if (ntokens != ncols) ret = 1;
      
      /*Fill row. */
      for (unsigned i=0; i<ncols; i++) {
	CSVTableColumn<Tvar>& col = columns[i];
	col.addRows(1);
	if (i < ntokens) {
	  std::string& v = tokens[i]; 
	  col[nrows] = v;
	}
      }
      nrows++;
    }
      
    return ret;
  }  

 protected:

  void split( std::vector<std::string> & theStringVector,
	      const  std::string  & theString,
	      const  std::string  & theDelimiter)
  {
    size_t  start = 0, end = 0;
    while ( end != std::string::npos)
    {
        end = theString.find( theDelimiter, start);
        // If at end, use length=maxLength.  Else use length=end-start.
        theStringVector.push_back( theString.substr( start,
                       (end == std::string::npos) ? std::string::npos : end - start));
        // If at end, use start=maxSize.  Else use start=end+delimiter.
        start = (   ( end > (std::string::npos - theDelimiter.size()) )
                  ?  std::string::npos  :  end + theDelimiter.size());
    }
  }

  // trim from both ends
  static inline std::string &trim(std::string &s) {
    
    /* If first is space, eat the consecutive spaces. */
    std::string::iterator first = s.begin();
    while((first != s.end()) && isspace(*first)) first++;
    if (first != s.begin()) s.erase(s.begin(), first-1);

    if (s.size() == 0) return s;

    /* If we reach this point, there is at least one non-space character. */
    /* If last is space, eat the consecutive spaces. */
    std::string::iterator last = s.end()-1;
    /* If first is space, eat the consecutive spaces. */
    while(isspace(*last)) last--;
    if (last != s.end()) s.erase(++last, s.end());

    return s;
  }

  /** Array of columns. */
  std::vector<CSVTableColumn<Tvar> > columns;
  
  unsigned nrows;
};

/** CSV String table. */ 
typedef CSVTable<std::string> CSVStringTable;

#endif // CSVTABLE_H

