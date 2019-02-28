// Name:     Disjunction.cpp
// Author:   A. M. Kazachkov
// Date:     2018-02-22
//-----------------------------------------------------------------------------
#include "Disjunction.hpp"
#include <limits>
#include <cmath> // abs

// Project files
#include "CglVPC.hpp" // ExitReason

/****************** PUBLIC  **********************/
/** Default constructor */
Disjunction::Disjunction() {
  initialize();
} /* default constructor */

/** Copy constructor */
Disjunction::Disjunction(const Disjunction& source) {
  initialize(&source);
} /* copy constructor */

/** Destructor */
Disjunction::~Disjunction() {
  for (DisjunctiveTerm term : terms) {
    if (term.basis) {
      delete term.basis;
      term.basis = NULL;
    }
  }
} /* destructor */

/** Assignment operator */
Disjunction& Disjunction::operator=(const Disjunction& source) {
  if (this != &source) {
    initialize(&source);
  }
  return *this;
} /* assignment operator */

/** Set up the disjunction class as new (except the timer pointer) */
void Disjunction::setupAsNew() {
  this->name = "";
  this->best_obj = std::numeric_limits<double>::max();
  this->worst_obj = std::numeric_limits<double>::lowest();
//  this->min_nb_obj_val = std::numeric_limits<double>::max();
  this->integer_obj = std::numeric_limits<double>::max();
  this->integer_sol.resize(0);
  this->common_changed_bound.resize(0);
  this->common_changed_value.resize(0);
  this->common_changed_var.resize(0);
  this->common_ineqs.resize(0);
  this->num_terms = 0;
  this->terms.resize(0);
} /* setupAsNew */

void Disjunction::setCgsName(std::string& cgsName,
    const std::string& disjTermName) {
  if (disjTermName.empty()) {
    return;
  }
  if (!cgsName.empty()) {
    cgsName += " V ";
  }
  cgsName += "(";
  cgsName += disjTermName;
  cgsName += ")";
} /* setCgsName (given disj term name) */

void Disjunction::setCgsName(std::string& cgsName, const int num_coeff,
    const int* const termIndices, const double* const termCoeff,
    const double termRHS, const bool append) {
if (num_coeff == 0) {
  return;
}
if (!cgsName.empty()) {
  if (!append) {
    cgsName += " V ";
  } else {
    cgsName.resize(cgsName.size() - 1); // remove last ")"
    cgsName += "; ";
  }
}
cgsName += append ? "" : "(";
for (int coeff_ind = 0; coeff_ind < num_coeff; coeff_ind++) {
  const double absCurrCoeff = std::abs(termCoeff[coeff_ind]);
  cgsName += (termCoeff[coeff_ind] > 0) ? "+" : "-";
  if (!isVal(absCurrCoeff, 1.)) {
    if (isVal(absCurrCoeff, std::floor(absCurrCoeff))
        || isVal(absCurrCoeff, std::ceil(absCurrCoeff)))
      cgsName += std::to_string(static_cast<int>(absCurrCoeff));
    else
      cgsName += std::to_string(absCurrCoeff);
  }
  cgsName += "x";
  cgsName += std::to_string(termIndices[coeff_ind]);
}
cgsName += " >= ";

if (isVal(termRHS, std::floor(termRHS))
    || isVal(termRHS, std::ceil(termRHS)))
  cgsName += std::to_string(static_cast<int>(termRHS));
else
  cgsName += std::to_string(termRHS);
cgsName += ")";
} /* setCgsName (one ineq per term) */

void Disjunction::setCgsName(std::string& cgsName, const int num_ineq_per_term,
    const std::vector<std::vector<int> >& termIndices,
    const std::vector<std::vector<double> >& termCoeff,
    const std::vector<double>& termRHS, const bool append) {
  if (num_ineq_per_term == 0) {
    return;
  }
  if (!cgsName.empty()) {
    if (!append) {
      cgsName += " V ";
    } else {
      cgsName.resize(cgsName.size() - 1); // remove last ")"
      cgsName += "; ";
    }
  }
  cgsName += append ? "" : "(";
  for (int i = 0; i < num_ineq_per_term; i++) {
    setCgsName(cgsName, termIndices[i].size(), termIndices[i].data(),
        termCoeff[i].data(), termRHS[i], (i > 0));
  }
  cgsName += ")";
} /* setCgsName */

/****************** PROTECTED **********************/
void Disjunction::initialize(const Disjunction* const source) {
  if (source != NULL) {
    this->name = source->name;
    this->best_obj = source->best_obj;
    this->worst_obj = source->worst_obj;
    this->integer_obj = source->integer_obj;
    this->integer_sol = source->integer_sol;
    this->timer = source->timer;
    this->common_changed_bound = source->common_changed_bound;
    this->common_changed_value = source->common_changed_value;
    this->common_changed_var = source->common_changed_var;
    this->common_ineqs = source->common_ineqs;
    this->num_terms = source->num_terms;
    this->terms = source->terms;
  } else {
    this->timer = NULL;
    setupAsNew();
  }
} /* initialize */

void Disjunction::updateObjValue(const double objVal) {
  if (objVal < this->best_obj) {
    this->best_obj = objVal;
  }
  if (objVal > this->worst_obj) {
    this->worst_obj = objVal;
  }
} /* updateObjValue */