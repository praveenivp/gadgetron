#pragma once
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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef MAX
    #define MAX(a,b)  ((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
    #define MIN(a,b)    ((a)>(b)?(b):(a))
#endif
#ifndef GRAD_RASTER_TIME
    #define GRAD_RASTER_TIME 10
#endif
namespace Gadgetron{

namespace Spiral{
	enum eSpiralType {
	SpiralOut = 1,
	SpiralIn = 2,
	DoubleSpiral = 3,
	SpiralInAndOut=4
};

class PESpiralDesign {
protected:
	std::vector<float> m_vfGx, m_vfGy, m_vfGz;
	double m_dAx, m_dAy, m_dAz, m_dAmp, m_dMomX, m_dMomY, m_dMomZ, m_dPreMomX, m_dPreMomY, m_dPreMomZ, m_dPostMomX, m_dPostMomY, m_dPostMomZ;
	double m_dGradRasterTime;
	double m_dLarmorConst;
	eSpiralType m_eSpiralType;

	double m_dResolution, m_dMaxAmplitude, m_dMinRiseTime;
	int m_Nitlv;
	std::vector<double> m_fov, m_radius;
	


public:

	// //constructer with ismrmrd header
	// PESpiralDesign(const ISMRMRD::AcquisitionHeader &acq_header)
	// {
		
	// }
	PESpiralDesign(void); // Constructor
	~PESpiralDesign(void); // Destructor

public:
	void setparameters(int Nitlv, double res, double* fov_, int nfov, double* radius_, int nradius, double Gmax, double Smax, eSpiralType spiralType, double T);
	void PrintOtherparameters();
	std::vector<float> getGx();
	std::vector<float> getGy();
	void printGX();
	void printGY();
	bool vdSpiralDesign(int Nitlv, double res, double* fov_, int nfov, double* radius_, int nradius, double Gmax, double Smax, eSpiralType spiralType, double T);
	bool setSpiralType(eSpiralType spiralType);
	bool calcTrajectory(std::vector<float> &vfKx, std::vector<float> &vfKy, std::vector<float> &vfDcf, long lADCSamples, int gridsize, double dADCshift, double dGradDelay);
	std::vector<float> jacksonDCF(std::vector<float> &vfKx, std::vector<float> &vfKy, int gridsize, float zeta);
};
}
}
