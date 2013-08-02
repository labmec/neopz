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
#include "arglib.h"

//#include <unordered_map>
//#define MAP unordered_map

#include <map>
#define MAP map

using namespace std;

//#define DEBUG_AL(...) cout << __VA_ARGS__
#define DEBUG_AL(...)

namespace clarg {

  /**
   * Arguments container. All the arguments are added into this object.
   */
  class args_container {
  public:

    void register_argument(arg_base* arg)
    {
      DEBUG_AL("Registering argument " << arg->get_name() << " (" 
	       << arg->get_desc() << ")" << endl);
      args[arg->get_name()] = arg;
    }

    void arguments_descriptions(ostream& os, string prefix, string suffix)
    {
      MAP<string, arg_base*>::iterator it;
      for (it = args.begin(); it != args.end(); it++) {
	arg_base* arg = it->second;
	os << prefix << arg->get_name() << " : " << arg->get_desc() << suffix;
      }
    }

    void list_arguments(ostream& os, bool defined_only)
    {
      MAP<string, arg_base*>::iterator it;
      for (it = args.begin(); it != args.end(); it++) {
	if (defined_only && !(it->second->was_set()))
	  continue;
	os << "# " << it->second->get_desc() << " (default value:";
	it->second->write_parameters(os, true);
	os << ")" << endl;
	os << it->second->get_name() << " ";
	it->second->write_parameters(os);
	os << endl;
      }
    }

    /**
    * Read the arguments from file.
    * Returns 0 if ok, != 0 otherwise.
    * TODO: this is not implemented yet. The goal of this method is to parser
    *       all the arguments from an argument file. In this way, you may dump
    *       all the arguments used in a given run and re-use them by parsing
    *       from the file. To do so, you may infoke this method.
    */
    int parse_arguments_from_file(std::istream& is)
    {
      /* Suggestion: parse lines, one by one. For each line:
       * - discard characters after '#'
       * - \# is converted into '#'
       * - build argc and argv for the line.
       * - call parse_arguments(argc, argv)
       */
      cerr << "ERROR: " << __PRETTY_FUNCTION__ << " was not implemented yet" << endl;
      return 1; /* Return error for now. */
    }
    
    int parse_arguments(int argc, char *argv[])
    {
      if (argc <= 0)
	return -1;
      prog_name = string(argv[0]);
      argv++;
      argc--;

      while (argc > 0) {
	
	string argn = string(argv[0]);
	/* Search argument */
	DEBUG_AL("parse_arguments: look for \"" << argn << "\"" << endl);
	std::MAP<string,arg_base*>::const_iterator got = args.find(argn);
	if (got == args.end()) {
	  cerr << "Error: could not find argument: " << argn << endl;
	  return -1;
	}
	DEBUG_AL("found it" << endl);
	arg_base* argp = got->second;
	argv++;
	argc--;

	/* Sanity checking. */
	if (argp->was_set()) {
	  /* argument is being set more than once. */
	  cerr << "WARNING: argument " << argp->get_name() << " is being set more than once." << endl;
	}
	else {
	  argp->mark_set(true);
	}

	/* Parse parameters */
	int consumed = argp->parse_parameters(argc, argv);
	DEBUG_AL("parse_arguments: consumed " << consumed 
		 << " parameters when parsing argument "
		 << argn << endl);
		 
	if (consumed < 0) {
	  /* parse_arguments returned error. */
	  cerr << "Error when parsing parameters for argument: " << got->first << endl;
	  return 1;
	}

	/* Update argc, argv */
	argv += consumed;
	argc -= consumed;
      }
      
      return 0;
    }

    /* Arguments. */
    MAP<string, arg_base*> args;
    /* Program name. */
    string prog_name;
  };

  class container_manager {
  public:
    container_manager() {create_container();}
    ~container_manager() {if (container != NULL) delete container; }

    args_container* get_container() { 
      if (container == NULL) 
	create_container();
      return container; 
    }

  private:

    static args_container* container;

    void create_container() {
      if (container == NULL)
	container = new args_container();
    }
  };

  args_container* container_manager::container = NULL;

  container_manager cm;


  arg_base::arg_base(const char* name, const char* desc) :
    arg_name (name), arg_desc (desc)
  {
    cm.get_container()->register_argument(this);
  }

  /* Print the arguments. */
  void values(ostream& os, bool defined_only)
  {
    cm.get_container()->list_arguments(os, defined_only);
  }
  
  /* Print the arguments descriptions. */
  void arguments_descriptions(ostream& os, string prefix, string suffix)
  {
    cm.get_container()->arguments_descriptions(os,prefix,suffix);
  }
  
  /* Parse the arguments. */
  int parse_arguments(int argc, char *argv[])
  {
    return cm.get_container()->parse_arguments(argc, argv);
  }

  /**
   * Dump the arguments to file.
   * Returns 0 if ok, != 0 otherwise.
   */
  int dump_arguments_to_file(ostream& of)
  {
    try {cm.get_container()->list_arguments(of, false);}
    catch (...){ return 1;}
    return 0;
  }

  /**
   * Read the arguments from file.
   * Returns 0 if ok, != 0 otherwise.
   */
  int parse_arguments_from_file(std::istream& is)
  {
    return cm.get_container()->parse_arguments_from_file(is);
  }

}

