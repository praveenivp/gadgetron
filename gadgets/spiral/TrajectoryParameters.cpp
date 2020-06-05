//
// Created by dchansen on 10/2/18.
//

#include <cpu/hoNDArray_utils.h>
#include <cpu/hoNDArray_fileio.h>
#include "TrajectoryParameters.h"
#include "vector_td_utilities.h"

namespace Gadgetron {
    namespace Spiral { 
        TrajectoryParameters::TrajectoryParameters(const ISMRMRD::IsmrmrdHeader &h) {
            ISMRMRD::TrajectoryDescription traj_desc;

            if (h.encoding[0].trajectoryDescription) {
                traj_desc = *h.encoding[0].trajectoryDescription;
            } else {
                throw std::runtime_error("Trajectory description missing");
            }

            if (traj_desc.identifier != "PESpiral") {
                throw std::runtime_error("Expected trajectory description identifier 'PESpiral', not found.");
            }


            try {
                auto userparam_long = to_map(traj_desc.userParameterLong);
                auto userparam_double = to_map(traj_desc.userParameterDouble);
                //piv edit
                Spiral_type=userparam_long.at("SpiralType");
                Resolution_mm=userparam_double.at("Resolution_mm");
                GradDelay_us=userparam_double.at("GradDelay_us");
                ADCShift_us=userparam_double.at("ADCShift_us");
                Tsamp_ns_ = userparam_long.at("DwellTime_ns");
                Nints_ = userparam_long.at("interleaves");


                gmax_ = userparam_double.at("MaxGradient_mT_per_m");
                smax_ = userparam_double.at("Slewmax_mT_m_ms");
                krmax_ = 1;//luserparam_double.at("krmax_per_cm");
                fov_ = h.encoding[0].reconSpace.fieldOfView_mm.x; //userparam_double.at("FOVCoeff_1_cm");
            } catch (std::out_of_range exception) {
                std::string s = "Missing user parameters: " + std::string(exception.what());
                throw std::runtime_error(s);

            }

            long ADCsamples=h.encoding[0].encodedSpace.matrixSize.x;
            TE_ = h.sequenceParameters->TE->at(0);


            if (h.userParameters) {
                try {
                    auto user_params_string = to_map(h.userParameters->userParameterString);
                    auto user_params_double = to_map(h.userParameters->userParameterDouble);

                    auto girf_kernel_string = user_params_string.at("GIRF_kernel");
                    this->girf_kernel = std::make_optional<hoNDArray<std::complex<float>>>(
                            GIRF::load_girf_kernel(girf_kernel_string));

                    girf_sampling_time_us = user_params_double.at("GIRF_sampling_time_us");

                } catch (std::out_of_range exception) { }
            }

            GDEBUG("GradDelay:                    %f\n", GradDelay_us);
            GDEBUG("ADCshift:                    %f\n", ADCShift_us);
            GDEBUG("ADCLength:                  %d\n", ADCsamples);
            GDEBUG("smax:                    %f\n", smax_);
            GDEBUG("gmax:                    %f\n", gmax_);
            GDEBUG("Tsamp_ns:                %d\n", Tsamp_ns_);
            GDEBUG("Nints:                   %d\n", Nints_);
            GDEBUG("fov:                     %f\n", fov_);
            GDEBUG("krmax:                   %f\n", krmax_);
            GDEBUG("Resolution_mm:                   %f\n",Resolution_mm);
            GDEBUG("SpiralType:                   %d\n", Spiral_type);
            GDEBUG("GIRF kernel:             %d\n", bool(this->girf_kernel));

        }
    

                std::pair<hoNDArray<floatd2>, hoNDArray<float>>
        TrajectoryParameters::calculate_trajectories_and_weight(const ISMRMRD::AcquisitionHeader &acq_header) {            	


            int nfov = 4;
            //peNC_spiral1_fa30_TR100_23052019.edb
            double* fov = new double[nfov];
            double* radius = new double[nfov];
            fov[0] = fov_; fov[1] = fov_; fov[2] = fov_; fov[3] = fov_;

            radius[0] = 0; radius[1] = 0.15; radius[2] = 0.2; radius[3] = 1;
           
            double dGradMaxAmpl = 40;
        
            Gadgetron::Spiral::eSpiralType typ = eSpiralType(Spiral_type);

            Gadgetron::Spiral::PESpiralDesign spiralTraj;
            
            Gadgetron::Spiral::PESpiralDesign  spiralTraj2;
            
            


            if (typ == SpiralInAndOut)
            {
                // SpiralInAndOut is not convertible to vdspiral::eSpiralType!
                spiralTraj.setparameters(Nints_, Resolution_mm, fov, nfov, radius, nfov,gmax_, smax_, eSpiralType::SpiralIn, GRAD_RASTER_TIME);
                spiralTraj.vdSpiralDesign(Nints_, Resolution_mm, fov, nfov, radius, nfov, gmax_, smax_, eSpiralType::SpiralIn, GRAD_RASTER_TIME);
            }
            else {
                spiralTraj.setparameters(Nints_, Resolution_mm, fov, nfov, radius, nfov, dGradMaxAmpl, smax_, typ, GRAD_RASTER_TIME);
                spiralTraj.vdSpiralDesign(Nints_, Resolution_mm, fov, nfov, radius, nfov, dGradMaxAmpl, smax_, typ, GRAD_RASTER_TIME);
            }

            

            if (typ == SpiralInAndOut) {
                // SpiralInAndOut is not convertible to vdspiral::eSpiralType!
                // spiral out
                spiralTraj.setSpiralType(eSpiralType::SpiralIn);

                // spiral in
                spiralTraj2 = spiralTraj;
                spiralTraj.setSpiralType(eSpiralType::SpiralOut);
            }
            else {
                spiralTraj.setSpiralType(typ);
            }
             
            long ADCsamples=(long)acq_header.number_of_samples;
            std::vector<float> kx;
            std::vector<float> ky;
            std::vector<float> dcf;
            spiralTraj.calcTrajectory(kx,ky,dcf,ADCsamples,std::ceil(fov_/Resolution_mm),ADCShift_us,GradDelay_us);


            hoNDArray<floatd2> trajectories((Nints_*ADCsamples));
            size_t nsamples = trajectories.get_number_of_elements();
             GDEBUG("nsamples:             %d\n", nsamples);
             GDEBUG("sizeof(dcf)/sizeof(dcf[0]) %f\n", (sizeof(dcf)/sizeof(dcf[0])));
            std::transform(kx.begin(),kx.begin()+nsamples,ky.begin(),trajectories.begin(),[](auto x, auto y){return floatd2(x,y)/10;});
            trajectories.reshape({ADCsamples,Nints_});
            hoNDArray<float> weights(nsamples);
            std::transform(dcf.begin(),dcf.begin()+nsamples,weights.begin(),[](auto x){return x;});
            weights.reshape({ADCsamples,Nints_});

            std::ofstream myfile1;
            myfile1.open ("/home/praveen/dcf.txt");
            for(std::vector<float>::iterator it = dcf.begin();it!=dcf.end();++it){
                myfile1 << *it<< ", "; //<<ky[i]<<", "<<dcf[i]<< ", ";
            }
            myfile1<<std::endl;
            myfile1.close();

            std::ofstream myfile2;
            myfile2.open ("/home/praveen/ky.txt");
            for(std::vector<float>::iterator it = ky.begin();it!=ky.end();++it){
                myfile2 << *it<< ", "; //<<ky[i]<<", "<<dcf[i]<< ", ";
            }
            myfile2<<std::endl;
            myfile2.close();
            std::ofstream myfile;
            myfile.open ("/home/praveen/kx.txt");
            for(std::vector<float>::iterator it = kx.begin();it!=kx.end();++it){
                myfile << *it<< ", "; //<<ky[i]<<", "<<dcf[i]<< ", ";
            }
            myfile<<std::endl;
            myfile.close();


            // if (this->girf_kernel){
            //     base_gradients = correct_gradients(base_gradients,Tsamp_ns_*1e-3,this->girf_sampling_time_us,acq_header.read_dir,acq_header.phase_dir,acq_header.slice_dir);
            //     //Weights should be calculated without GIRF corrections according to Hoge et al 2005
            //     trajectories = calculate_trajectories(base_gradients,sample_time,krmax_);
            // }

            return std::make_pair(std::move(trajectories), std::move(weights));

        }

        hoNDArray<floatd2>
        TrajectoryParameters::correct_gradients(const hoNDArray<floatd2> &gradients, float grad_samp_us,
                                                float girf_samp_us, const float *read_dir, const float *phase_dir,
                                                const float *slice_dir) {

            arma::fmat33 rotation_matrix;
            rotation_matrix(0, 0) = read_dir[0];
            rotation_matrix(0, 1) = read_dir[1];
            rotation_matrix(0, 2) = read_dir[2];
            rotation_matrix(1, 0) = phase_dir[0];
            rotation_matrix(1, 1) = phase_dir[1];
            rotation_matrix(1, 2) = phase_dir[2];
            rotation_matrix(2, 0) = slice_dir[0];
            rotation_matrix(2, 1) = slice_dir[1];
            rotation_matrix(2, 2) = slice_dir[2];

            return GIRF::girf_correct(gradients, *girf_kernel, rotation_matrix, grad_samp_us, girf_samp_us, TE_);

        }
    }
}