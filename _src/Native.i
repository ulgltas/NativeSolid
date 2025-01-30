// SWIG input file of the 'NativeSolid' module

%feature("autodoc","1");

%module(docstring=
"'NativeSolid' module",
directors="1",
threads="1"
) NativeSolid
%{

#include <string>
#include "NativeSolidSolver.h"

%}

%include "std_string.i"
%include "NativeSolidSolver.h"
