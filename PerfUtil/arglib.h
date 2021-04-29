/***************************************************************************
 *   Copyright (C) 2012 by Edson Borin                                     *
 *   edson@ic.unicamp.br                                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/**
 * This is a C++ command line argument library. The goal is to enable
 * the user to easily read arguments from the command line. The library
 * works as follows.
 *
 * 1) include the "arglib.h" file
 * 2) add a call to clarg::parse_arguments(arc, argv);
 * 3) create new arguments by instantiating objects from the clarg::argXXX classes.
 *    IMPORTANT: the objects must be declared on the global scope so that thei are
 *    created before the main function starts executing.
 *
 * EXAMPLE:
 *
 * #include "arglib.h"
 * 
 * clarg::argInt    n1("-n1", "integer argument", 0);
 * clarg::argDouble d1("-d1", "double argument", 0);
 * clarg::argString s1("-s1", "string argument", "default value");
 * clarg::argBool help("-h",  "display the help message");
 * 
 * int main(int argc, char *argv[])
 * {
 *   // Parse the arguments
 *   if (clarg::parse_arguments(argc, argv)) {
 *     cerr << "Error when parsing the arguments!" << endl;
 *     return 1;
 *   }
 * 
 *   std::cout << "Values set by user:"<< std::endl;
 *   clarg::values(cout);
 * 
 *   if (help.get_value() == true) {
 *     cout << "Usage: " << argv[0] << endl;
 *     clarg::arguments_descriptions(cout, "  ", "\n");
 *     return 1;
 *   }
 * 
 *   cout << "The user has provided the following arguments:" << endl;
 *   if (s1.was_set()) {
 *     cout << s1.get_name() << " = " << s1.get_value() << endl;
 *   }
 *   if (d1.was_set()) {
 *     cout << d1.get_name() << " = " << d1.get_value() << endl;
 *   }
 *   if (n1.was_set()) {
 *     cout << n1.get_name() << " = " << n1.get_value() << endl;
 *   }
 * 
 *   return 0;
 * }
 *
 * HINT: you may declare new arguments in any file, as long 
 *       as you call the clarg::parse_arguments function at the 
 *       beginning of the main function.
 *
 * NOTE: Notice that it is up to the user to verify if the required 
 *       arguments were provided. You may use the "was_set()" method.
 */
#ifndef ARGLIB_H
#define ARGLIB_H

#include <string>
#include <sstream>
#include <iostream>
#include <vector>

#include "pzreal.h" // PRETTY_FUNCTION for WIN32 Macro definition

namespace clarg 
{

  /** 
   * Parse the command line arguments. This must be called at the main.
   * Returns 0 if ok, != 0 otherwise.
   */
  int parse_arguments(int argc, char *argv[]);

  /**
   * Dump the arguments to file.
   * Returns 0 if ok, != 0 otherwise.
   */
  int dump_arguments_to_file(std::ostream& os);

  /**
   * Read the arguments from file.
   * Returns 0 if ok, != 0 otherwise.
   */
  int parse_arguments_from_file(std::istream& is);

  /**
   * Print the arguments descriptions. Very usefull when printing usage
   * messages.
   */
  void arguments_descriptions(std::ostream& os, std::string prefix, std::string suffix);

  /** 
   * Print the arguments values. If defined_only is true, only the arguments set
   * at the command line are printed.
   */
  void values(std::ostream& os, bool defined_only=true);

  class args_container;

  /** The base argument class. Manages arguments names and descriptions. */
  class arg_base
  {
  public:

    arg_base(const char* name, const char* desc);

    const std::string& get_name() const { return arg_name; }

    const std::string& get_desc() const { return arg_desc; }

    /* Return true if the argument was set, false otherwise. */
    bool was_set() const { return arg_set; }

    void mark_set(bool m) { arg_set = m; }

  protected:

    /** 
     * Parse the argument parameters. Must be implemented by the specialized
     * argument class.  
     * When parsing, make sure you copy the parsed parameters into the arg_parameters array.
     * Returns the number of parameters parsed if ok, -1 if an error occured.
     */
    virtual int parse_parameters (int argc, char* argv []) = 0;

    /**
     * write the argument parameters into the output stream.
     */
    virtual void write_parameters (std::ostream& os, bool def=false) const = 0;

    /** True if the argument was set */
    bool   arg_set;
    /** The argument name. */
    std::string arg_name;
    /** The argument description. */
    std::string arg_desc;

    friend class args_container;
  };

  /** The argument template class. */ 
  template<class T>
    class argT : public arg_base
  {
  public:

  argT(const char* arg, const char* desc) :
    arg_base (arg, desc) {}

    /* Returns the argument value. */ 
    const T& get_value() const { return value; }

    /* Sets the argument value. */ 
    void set_value(const T& v) const { value = v; }

  protected:

    /* Argument value */
    T value;
    /* Default argument value */
    T def_value;
  };

  /** 
   * String argument class.
   */
  class argString : public argT<std::string>
  {
  public:
  argString(const char* arg, const char* desc, std::string v) :
    argT<std::string>(arg,desc)
  {
    def_value = v;
    value = def_value;
  }
  protected:
    int parse_parameters (int argc, char* argv []) override
    {
      if (argc <= 0) return -1;
      value = std::string(argv[0]);
      return 1;
    }
    void write_parameters (std::ostream& os, bool def) const override
    {
      if (def)
	os << def_value;
      else
	os << value;
    }

  };

  /** 
   * Integer argument class.
   */
  class argInt : public argT<int>
  {
  public:
  argInt(const char* arg, const char* desc, int v = 0) :
      argT<int>(arg,desc)
    {
      def_value = v;
      value = def_value;
    }
  protected:
    int parse_parameters (int argc, char* argv []) override
    {
      if (argc <= 0) return -1;
      try {
	std::stringstream(argv[0]) >> value;
      }
      catch (const std::exception& ) {return -1;}
      return 1;
    }
    void write_parameters (std::ostream& os, bool def) const override
    {
      if (def)
	os << def_value;
      else
	os << value;
    }
  };

  /** 
   * Double argument class.
   */
  class argDouble : public argT<double>
  {
  public:
  argDouble(const char* arg, const char* desc, double v = 0.0) :
      argT<double>(arg,desc)
    {
      def_value = v;
      value = def_value;
    }
  protected:
    int parse_parameters (int argc, char* argv []) override
    {
      if (argc <= 0) return -1;
      try {
	std::stringstream(argv[0]) >> value;
      }
      catch (const std::exception& ) {return -1;}
      return 1;
    }
    void write_parameters (std::ostream& os, bool def) const override
    {
      if (def)
	os << def_value;
      else
	os << value;
    }
  };

  /** 
   * Boolean argument class.
   */
  class argBool : public argT<bool>
  {
  public:
    argBool(const char* arg, const char* desc, bool v = false) :
      argT<bool>(arg,desc)
    {
      def_value = v;
      value = def_value;
    }
  protected:
    int parse_parameters (int argc, char* argv []) override
    { value = true; return 0; }
    void write_parameters (std::ostream& os, bool def) const override
    {}
  };

}

#endif
