#include <iostream>
//#include "D:\VM Shared\_cmath.h"
//#include "D:\VM Shared\_vector.h"
#include "spline.h"
#include "mtg_functions_pe.h"
#include<vector>

#include <ismrmrd/xml.h>
#include "log.h"
#include "Gadget.h"
#include "mri_core_utility.h"


//only for gnuplot
// Warn about use of deprecated functions.
#define GNUPLOT_DEPRECATE_WARN
//#include "gnuplot-iostream.h"
#include <map>
#include <limits>
#include <cmath>
#include<fstream>


enum eSpiralType {
	SpiralOut = 1,
	SpiralIn = 2,
	DoubleSpiral = 3,
	SpiralInAndOut=4
};

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

