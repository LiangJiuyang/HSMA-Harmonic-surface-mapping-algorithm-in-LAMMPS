/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_ERROR_H
#define LMP_ERROR_H

#include "pointers.h"

#ifdef LAMMPS_EXCEPTIONS
#include "exceptions.h"
#endif

namespace LAMMPS_NS {

class Error : protected Pointers {
 public:
  Error(class LAMMPS *);

  [[ noreturn ]] void universe_all(const std::string &, int, const std::string &);
  [[ noreturn ]] void universe_one(const std::string &, int, const std::string &);
  void universe_warn(const std::string &, int, const std::string &);

  [[ noreturn ]] void all(const std::string &, int, const std::string &);
  [[ noreturn ]] void one(const std::string &, int, const std::string &);
  void warning(const std::string &, int, const std::string &, int = 1);
  void message(const std::string &, int, const std::string &, int = 1);
  [[ noreturn ]] void done(int = 0); // 1 would be fully backwards compatible

#ifdef LAMMPS_EXCEPTIONS
  std::string get_last_error() const;
  ErrorType get_last_error_type() const;
  void set_last_error(const std::string &msg, ErrorType type = ERROR_NORMAL);

 private:
  std::string last_error_message;
  ErrorType last_error_type;
#endif
};

}

#endif

/* ERROR/WARNING messages:

*/
