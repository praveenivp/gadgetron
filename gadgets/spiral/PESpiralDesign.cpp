#include "PESpiralDesign.h"

namespace Gadgetron
{
namespace Spiral
{

	// //constructer with ismrmrd header
	// PESpiralDesign(const ISMRMRD::AcquisitionHeader &acq_header)
	// {
		
	// }
	PESpiralDesign::PESpiralDesign(void) // Constructor
	{
		this->m_eSpiralType = SpiralOut;
		this->m_Nitlv = 8;
		this->m_dResolution = 2.;//mm
		this->m_dLarmorConst = 42.575575;//42.5756;
	}


	PESpiralDesign::~PESpiralDesign(void) // Destructor
	{}


	void PESpiralDesign::setparameters(int Nitlv, double res, double* fov_, int nfov, double* radius_, int nradius, double Gmax, double Smax, eSpiralType spiralType, double T)
	{
		m_Nitlv = Nitlv;
		m_dResolution = res;
		m_dMaxAmplitude = Gmax;
			m_dMinRiseTime=  1000. / Smax;
			m_eSpiralType = spiralType;
			m_dGradRasterTime = T;
			m_fov.clear();
			m_radius.clear();
			for (int i = 0; i < nfov; i++)
			{
				m_fov.push_back(fov_[i]);
				m_radius.push_back(radius_[i]);
			}
			m_dLarmorConst = 42.575575;

	}

	void PESpiralDesign::PrintOtherparameters()
	{
		std::cout << this->m_dAx << std::endl;
		std::cout << this->m_dMomX << std::endl;
		std::cout << this->m_dMomY << std::endl;
		std::cout << this->m_dMomZ << std::endl;
		std::cout << this->m_dPreMomX << std::endl;
		std::cout << this->m_dPreMomY << std::endl;
		std::cout << this->m_dPreMomZ << std::endl;
		std::cout << this->m_dPostMomX << std::endl;
		std::cout << this->m_dPostMomY << std::endl;
		std::cout << this->m_dPostMomZ << std::endl;
		

	}

	std::vector<float> PESpiralDesign::getGx()
	{
		return m_vfGx;
	}
	std::vector<float> PESpiralDesign::getGy()
	{
		return m_vfGy;
	}
	void PESpiralDesign::printGX()
	{
		
		for (int i=0;i<this->m_vfGx.size();i++)
			std::cout << m_vfGx[i]<< "  ";
		std::ofstream dumpfile;
		std::string filename = "C:\\Users\\pvalsala\\Documents\\MATLAB\\dumpfile.txt";
		dumpfile.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);
			for (int i = 0; i < this->m_vfGx.size(); i++)
				dumpfile << m_vfGx[i] *m_dAx<< "     " <<m_vfGy[i]*m_dAy<<std::endl;
			dumpfile.close();
			
		
	}

	void PESpiralDesign::printGY()
	{
		for (int i = 0; i < this->m_vfGy.size(); i++)
			std::cout << m_vfGy[i] << "  ";
	}
	
	bool PESpiralDesign::vdSpiralDesign(int Nitlv, double res, double* fov_, int nfov, double* radius_, int nradius, double Gmax, double Smax, eSpiralType spiralType, double T) {

		long k; // loop index

		if (nfov != nradius) {
			std::cout << "Error: array fov needs to have same size as array radius" << std::endl;
			return false;
		}
		// we need to make copies of fov and radius since we do a unit transformation
		double *fov = new double[nfov];
		double *radius = new double[nradius];

		double kmax = 5. / res;  // kmax = 1/(2*res) BUT: kmax in 1/cm, res in mm

		double fovmax = 0.;
		for (k = 0; k < nfov; ++k) {
			fov[k] = fov_[k] / 10.; // mm->cm
			if (fov[k] > fovmax)
				fovmax = fov[k];
			radius[k] = kmax * radius_[k];
		}


		//double dr = 1./500. * 1./(fovmax/Nitlv);
		double dr = 1. / 100. * 1. / (fovmax / Nitlv); // a little faster
		long   nr = long(kmax / dr) + 1;

		std::vector<double> x, y, z;
		x.resize(nr, 0.);
		y.resize(nr, 0.);
		z.resize(nr, 0.);

		double theta = 0.;
		for (k = 0; k < nr; k++) {
			double r = k * dr;
			double cFoV = fov[nfov - 1];
			for (int l = 0; l < nfov; ++l) {
				if (r < radius[l]) {
					if (l == 0 || l == nfov - 1)
						cFoV = fov[l];
					else {// linearer übergang // linear transition
						double step = (r - radius[l - 1]) / (radius[l] - radius[l - 1]);
						cFoV = step * fov[l] + (1. - step)*fov[l - 1];
					}
					break;
				}
			}
			x[k] = r * cos(theta);
			y[k] = r * sin(theta);
			if (spiralType == DoubleSpiral)
				theta += M_PI * dr*cFoV / Nitlv;
			else
				theta += 2.*M_PI*dr*cFoV / Nitlv;
		}

		delete[] fov; delete[] radius;

		Gmax /= 10.;   // mT/m    -> G/cm
		Smax /= 10.;   // mT/m/ms -> G/cm/ms
		T /= 1000.; // us      -> ms

		int n;
		double g0 = 0.; // to simplify sequence development, our gradient will start at 0.
		double gfin = 0.; // and end at 0.
		double *gx; double *gy; double *gz;
		//clock_t start = clock();
		minTimeGradientRIV(&x[0], &y[0], &z[0], nr, g0, gfin, Gmax, Smax, T, gx, gy, gz, n, -1.5, m_dLarmorConst / 10.);
		/*gx = &x[0];
		gy = &y[0];
		gz = &z[0];
		n = nr;*/
		//clock_t end = clock();
		//cout << "start = " << start << "  end = " << end << "  end-start = " << end-start << "  CLOCKS_PER_SEC = " << CLOCKS_PER_SEC << endl;

		// determine max gradient amplitudes
		m_dAx = 0.;
		m_dAy = 0.;
		m_dAmp = 0.;
		for (k = 0; k < n; ++k) {
			if (fabs(gx[k]) > m_dAx)
				m_dAx = fabs(gx[k]);
			if (fabs(gy[k]) > m_dAy)
				m_dAy = fabs(gy[k]);
			double dAmp = sqrt(gx[k] * gx[k] + gy[k] * gy[k]);
			if (dAmp > m_dAmp)
				m_dAmp = dAmp;
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////
		///// pehses: start gradient with zero and ramp down the gradient to zero for good at the end /////
		///////////////////////////////////////////////////////////////////////////////////////////////////
		int addRampUp = 1;
		int addRampDn = 1 + (int)(sqrt(gx[n - 1] * gx[n - 1] + gy[n - 1] * gy[n - 1] + gz[n - 1] * gz[n - 1]) / Smax / T);
		if (spiralType == DoubleSpiral) {
			addRampUp = addRampDn;
			m_vfGx.resize(2 * n + addRampUp + addRampDn + 1, 0.);
			m_vfGy.resize(2 * n + addRampUp + addRampDn + 1, 0.);
			// first value of gradients is zero
			// copy & scale gradient to interval -1...+1
			for (k = 0; k < n; ++k) {
				m_vfGx[k + addRampUp] = (float)(gx[n - k - 1] / (m_dAx > 0 ? m_dAx : 1.));
				m_vfGy[k + addRampUp] = (float)(gy[n - k - 1] / (m_dAy > 0 ? m_dAy : 1.));
				m_vfGx[n + k + addRampUp + 1] = (float)(gx[k] / (m_dAx > 0 ? m_dAx : 1.));
				m_vfGy[n + k + addRampUp + 1] = (float)(gy[k] / (m_dAy > 0 ? m_dAy : 1.));
			}
			n = 2 * n + 1;
		}
		else {
			m_vfGx.resize(n + addRampUp + addRampDn, 0.);
			m_vfGy.resize(n + addRampUp + addRampDn, 0.);
			// first value of gradients is zero
			// copy & scale gradient to interval -1...+1
			for (k = 0; k < n; ++k) {
				m_vfGx[k + addRampUp] = (float)(gx[k] / (m_dAx > 0 ? m_dAx : 1.));
				m_vfGy[k + addRampUp] = (float)(gy[k] / (m_dAy > 0 ? m_dAy : 1.));
			}
		}
		delete[] gx; delete[] gy; delete[] gz;

		// linear ramp up from zero (we should already be close)
		for (k = 0; k < addRampUp; ++k) {
			m_vfGx[k] = (float)(m_vfGx[addRampUp] * (double(k) / addRampUp));
			m_vfGy[k] = (float)(m_vfGy[addRampUp] * (double(k) / addRampUp));
		}
		// linear ramp down to zero at the end (we should already be close)
		for (k = 0; k < addRampDn; ++k) {
			m_vfGx[n + addRampUp + k] = (float)(m_vfGx[n + addRampUp - 1] * (1. - (k + 1.) / addRampDn));
			m_vfGy[n + addRampUp + k] = (float)(m_vfGy[n + addRampUp - 1] * (1. - (k + 1.) / addRampDn));
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////

		// G/cm -> mT/m   
		m_dAx *= 10.;
		m_dAy *= 10.;
		m_dAmp *= 10.;

		// now calculate the gradient moments
		m_dMomX = 0.; m_dMomY = 0.; m_dMomZ = 0.;
		for (k = 0; k < (int)m_vfGx.size(); ++k) {
			m_dMomX += m_dAx * m_vfGx[k] * m_dGradRasterTime;
			m_dMomY += m_dAy * m_vfGy[k] * m_dGradRasterTime;
		}
		m_dPreMomX = 0.; m_dPreMomY = 0.; m_dPreMomZ = 0.;
		m_dPostMomX = 0.; m_dPostMomY = 0.; m_dPostMomZ = 0.;
		if (spiralType == SpiralIn) {
			m_dPreMomX = m_dMomX;
			m_dPreMomY = m_dMomY;
			m_dPreMomZ = m_dMomZ;
			// we have to time-reverse the trajectory!
			for (long k = 0; k < (int)m_vfGx.size() / 2; ++k) {
				std::swap(m_vfGx[k], m_vfGx[m_vfGx.size() - 1 - k]);
				std::swap(m_vfGy[k], m_vfGy[m_vfGy.size() - 1 - k]);
			}
		}
		else if (spiralType == SpiralOut) {
			m_dPostMomX = m_dMomX;
			m_dPostMomY = m_dMomY;
			m_dPostMomZ = m_dMomZ;
		}
		else if (spiralType == DoubleSpiral) {
			for (k = 0; k < (int)(m_vfGx.size() / 2); ++k) {
				m_dPreMomX += m_dAx * m_vfGx[k] * m_dGradRasterTime;
				m_dPreMomY += m_dAy * m_vfGy[k] * m_dGradRasterTime;
			}
			m_dPostMomX = m_dMomX - m_dPreMomX;
			m_dPostMomY = m_dMomY - m_dPreMomY;
			m_dPostMomZ = m_dMomZ - m_dPreMomZ;

		}
		else
			return false;

		return true;
	}


	bool PESpiralDesign::setSpiralType(eSpiralType spiralType) {
		bool bStatus = true;
		if (spiralType != m_eSpiralType) {
			if ((m_eSpiralType != DoubleSpiral) && (spiralType != DoubleSpiral)) {
				// we have to time-reverse the trajectory!
				for (long k = 0; k < (int)m_vfGx.size() / 2; ++k) {
					std::swap(m_vfGx[k], m_vfGx[m_vfGx.size() - 1 - k]);
					std::swap(m_vfGy[k], m_vfGy[m_vfGy.size() - 1 - k]);
				}
				std::swap(m_dPreMomX, m_dPostMomX);
				std::swap(m_dPreMomY, m_dPostMomY);
				std::swap(m_dPreMomZ, m_dPostMomZ);
				#ifdef BUILD_SEQU
				bStatus = vdspiral::prepGradients();
				#endif
				return bStatus;
			}
			else { // we need to recalculate the trajectory
				m_eSpiralType = spiralType;
				return false; //this->vdSpiralDesign(m_Nitlv, m_dResolution, m_fov, m_radius, m_dMaxAmplitude, m_dMinRiseTime, m_eSpiralType, m_dLarmorConst, m_dGradRasterTime);
			}
		}
		return true;
	}


bool PESpiralDesign::calcTrajectory(std::vector<float> &vfKx, std::vector<float> &vfKy, std::vector<float> &vfDcf, long lADCSamples, int gridsize, double dADCshift, double dGradDelay) {
    //only spiral out supported for now!
    
    if (m_vfGx.size()==0 || m_vfGx.size()!=m_vfGy.size())
        return false;
    
    long lGradSamples = m_vfGx.size();
    long k,l;
    
    double dwelltime = (lGradSamples*m_dGradRasterTime + dADCshift)/lADCSamples;
    // iflag used in spline method to signal error 
    int *iflag;
    int iflagp;
    iflag = &iflagp;
    int *last;
    int lastp = 0;
    last = &lastp;
    
    // kgradx,kgrady: k-space trajectory on gradient raster
    int nFillpre  = 2;
    int nFillpost = 2;
    if (m_eSpiralType == SpiralOut)
        nFillpre  += int(dADCshift/m_dGradRasterTime);
    else
        nFillpost += int(dADCshift/m_dGradRasterTime);
    
    long lFilledSamples = lGradSamples+nFillpre+nFillpost;
    double *kgradx  = new double[lFilledSamples];
    double *kgrady  = new double[lFilledSamples];
    for (k=0;k<nFillpre+1;++k) {
        kgradx[k] = 0.;
        kgrady[k] = 0.;
    }
    double cumsumx=0.,cumsumy=0.;
    for (k=1;k<lGradSamples;++k) {
        cumsumx += m_dAx * (m_vfGx[k]+m_vfGx[k-1])/2.;
        kgradx[k+nFillpre]  = cumsumx * m_dGradRasterTime * m_dLarmorConst/1e5;
        cumsumy += m_dAy * (m_vfGy[k]+m_vfGy[k-1])/2.;
        kgrady[k+nFillpre]  = cumsumy * m_dGradRasterTime * m_dLarmorConst/1e5;
    }
    for (k=0;k<nFillpost;++k) {
        cumsumx += m_dAx * m_vfGx[lGradSamples-1]/2.;
        cumsumy += m_dAy * m_vfGy[lGradSamples-1]/2.;
        kgradx[k+nFillpre+lGradSamples] = cumsumx * m_dGradRasterTime * m_dLarmorConst/1e5;
        kgrady[k+nFillpre+lGradSamples] = cumsumy * m_dGradRasterTime * m_dLarmorConst/1e5;
    }
    
    double *coeff1x = new double[lFilledSamples];
    double *coeff2x = new double[lFilledSamples];
    double *coeff3x = new double[lFilledSamples];
    double *coeff1y = new double[lFilledSamples];
    double *coeff2y = new double[lFilledSamples];
    double *coeff3y = new double[lFilledSamples];
    double *t_grad  = new double[lFilledSamples];
    // Initalize gradient raster time
    for (k=0;k<lFilledSamples;++k)
        t_grad[k] = (k-nFillpre+1)*m_dGradRasterTime;
    
    spline(lFilledSamples, 0, 0, 0, 0, t_grad, kgradx, coeff1x, coeff2x, coeff3x, iflag);
    spline(lFilledSamples, 0, 0, 0, 0, t_grad, kgrady, coeff1y, coeff2y, coeff3y, iflag);

    // --------------------------------------------------------------
    // Interpolated curve
    // --------------------------------------------------------------
    vfKx.resize(lADCSamples*m_Nitlv,0.);
    vfKy.resize(lADCSamples*m_Nitlv,0.);
    for (k=0; k<lADCSamples; k++) {
        // Time for ACD sampling point
        double t_relativeToGrad = (k+0.5)*dwelltime + dGradDelay;
        if (m_eSpiralType == SpiralOut)
            t_relativeToGrad -= dADCshift;
        else
            t_relativeToGrad += dADCshift;
        vfKx[k] = (float) seval(lFilledSamples, t_relativeToGrad, t_grad, kgradx, coeff1x, coeff2x, coeff3x, last);
        vfKy[k] = (float) seval(lFilledSamples, t_relativeToGrad, t_grad, kgrady, coeff1y, coeff2y, coeff3y, last);
    }
    delete[] kgradx;
    delete[] kgrady;
    delete[] t_grad;
    delete[] coeff1x;
    delete[] coeff2x;
    delete[] coeff3x;
    delete[] coeff1y;
    delete[] coeff2y;
    delete[] coeff3y;
    
    if (m_eSpiralType == SpiralIn) {
        //for spiral in: make sure that trajectory ends in the center of k-space
        for (k=0; k<lADCSamples; k++) {
            vfKx[k] -= vfKx[lADCSamples-1];
            vfKy[k] -= vfKy[lADCSamples-1];
        }
    } else if (m_eSpiralType == DoubleSpiral) {
        //for double spiral: make sure that middle of trajectory is in the center of k-space
        float midX = vfKx[lADCSamples/2];
        float midY = vfKy[lADCSamples/2];
        for (k=0; k<lADCSamples; k++) {
            vfKx[k] -= midX;
            vfKy[k] -= midY;
        }
    }
    
    //now calculate trajectory for other interleaves by rotation
    for (k=1;k<m_Nitlv;++k) {
        float phi = (float) (2.* M_PI * k / m_Nitlv);
        if (m_eSpiralType == DoubleSpiral) // we only need to distribute the spirals over M_PI 
            phi /= 2.f;
        for (l=0; l<lADCSamples; ++l) {
//             vfKx[l+k*lADCSamples] = (float) (cos(phi) * vfKx[l] + sin(phi) * vfKy[l]);
//             vfKy[l+k*lADCSamples] = (float) (-sin(phi) * vfKx[l] + cos(phi) * vfKy[l]);
            vfKx[l+k*lADCSamples] = (float) (cos(phi) * vfKx[l] - sin(phi) * vfKy[l]);
            vfKy[l+k*lADCSamples] = (float) (sin(phi) * vfKx[l] + cos(phi) * vfKy[l]);
        }
    }
    
    // Calculate density compensation function
    vfDcf = jacksonDCF(vfKx, vfKy, gridsize, 1.f);
    
    return true;
}


std::vector<float> PESpiralDesign::jacksonDCF(std::vector<float> &vfKx, std::vector<float> &vfKy, int gridsize, float zeta) {
    
    int k,l,m;
    long nsamples = vfKx.size()/m_Nitlv;
    
    //scale zeta:
    //find maxk:
    float kmax = 0.f;
    for (k=0;k<nsamples;++k) {
        float tmp = vfKx[k]*vfKx[k]+vfKy[k]*vfKy[k];
        if (tmp>kmax)
            kmax = fabs(vfKx[k]);
    }
    kmax  = sqrt(kmax);
    zeta *= 2.f * kmax * 2.f/gridsize;
    
    //cut dcf at 0.85 of kmax
    for (k=0;k<nsamples;++k) {
        float tmp = sqrt(vfKx[k]*vfKx[k]+vfKy[k]*vfKy[k]);
        if (tmp>0.85f*kmax)
            break;
    }
    // wi is cutoff and normalized to average wi between cutoff_ix1 and cutoff_ix2
    int cutoff_ix1 = k;
    int cutoff_ix2 = cutoff_ix1 + (nsamples-cutoff_ix1)/4+1; 
    
    float zeta_sq = zeta*zeta;
    std::vector<float> wi(nsamples*m_Nitlv,1.);
    for (k=0;k<nsamples;++k) {
        float goal = 0.;
        float kxk = vfKx[k];
        float kyk = vfKy[k];
        // vorsicht: Skript nimmt gleichverteilte Interleaves an (Winkel)
        for (l=0; l<m_Nitlv;l++) {
            for (m=0;m<nsamples;++m) {
                float dx = vfKx[m+l*nsamples] - kxk;
                dx = dx*dx;
                if (dx < zeta_sq) {
                    float dr = vfKy[m+l*nsamples] - kyk;
                    dr = dx + dr*dr;
                    if (dr < zeta_sq) {
                        dr = (float) sqrt(dr);
                        //simple hann filter
                        float kern = 0.5f - 0.5f * (float)cos(2.*M_PI*(1.+dr/zeta)/2.);
                        goal += kern;
                    }
                }
            }
        }
        wi[k] = 1.f/goal;
    }
    
    // determine cutoff value for wi
    float cutoff_val = 0.;
    for (k=cutoff_ix1;k<cutoff_ix2;++k) 
        cutoff_val += wi[k];
    cutoff_val /= MAX(1,cutoff_ix2-cutoff_ix1);
    // normalize wi by cutoff value
    for (k=0;k<nsamples;++k)
        wi[k] = MIN(1.f,wi[k]/cutoff_val);
            
    //now copy wi from interleave 0 to other interleaves;
    for (k=1; k<m_Nitlv;++k) {
        for (l=0; l<nsamples; ++l)
            wi[l+k*nsamples] = wi[l];
    }
    return wi;
}


}
}