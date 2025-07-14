#include "Observables_MultiOrbSF.h"

/*
 * ***********
 *  Functions in Class Observables ------
 *  ***********
*/


void Observables_MultiOrbSF::Get_red_den_mat(Mat_4_Complex_doub &Red_Den_mat, double mu){

    double Beta=Parameters_.beta;
    //Red_Den_mat; //c^dagger[j][orb+n_orb*spin]c[l][orb+n_orb*spin2]
    int c1, c2;
    int jx,jy, lx,ly;
    for (int j=0;j<Parameters_.ns;j++){
        jx = Coordinates_.indx_cellwise(j);
        jy = Coordinates_.indy_cellwise(j);

        for (int l=0;l<Parameters_.ns;l++){
            lx = Coordinates_.indx_cellwise(l);
            ly = Coordinates_.indy_cellwise(l);


            for (int orb=0;orb<Parameters_.n_orbs;orb++){
                for (int orb2=0;orb2<Parameters_.n_orbs;orb2++){

                    for (int s=0;s<2;s++){
                        for (int s2=0;s2<2;s2++){

                            Red_Den_mat[j][orb+Parameters_.n_orbs*s][l][orb2+Parameters_.n_orbs*s2]=zero_complex;
                            for (int n=0;n<Hamiltonian_.eigs_.size();n++){

                                c1 = Coordinates_.Nbasis_[jx][jy][orb] + s*(Coordinates_.nbasis_);
                                c2 = Coordinates_.Nbasis_[lx][ly][orb2] + + s2*(Coordinates_.nbasis_);

                                //remember psi goes with anhilation
                                //Red_Den_mat[j][s][l][s2] = <C_{l,s2}^{dagger} C_{j,s}>
                                Red_Den_mat[j][orb+Parameters_.n_orbs*s][l][orb2+Parameters_.n_orbs*s2] += conj(Hamiltonian_.Ham_(c2,n))*(Hamiltonian_.Ham_(c1,n))*(1.0/(1.0 + exp(Beta*(Hamiltonian_.eigs_[n] - mu))));

                            }

                        }
                    }

                }
            }
        }
    }


}


void Observables_MultiOrbSF::Get_Non_Interacting_dispersion()
{
}

double Observables_MultiOrbSF::Lorentzian(double x, double brd)
{
    double temp;

    temp = (1.0 / PI) * ((brd / 2.0) / ((x * x) + ((brd * brd) / 4.0)));

    return temp;
}

void Observables_MultiOrbSF::DensityOfStates()
{
    //-----------Calculate Bandwidth------//
    BandWidth = 2.0;
    //-----------------------------------//

} // ----------

void Observables_MultiOrbSF::OccDensity()
{

} // ----------

void Observables_MultiOrbSF::TotalOccDensity()
{

} // ----------

complex<double> Observables_MultiOrbSF::SiSjQ(int i, int j) { return SiSjQ_(i, j); }

double Observables_MultiOrbSF::SiSj(int i, int j) { return SiSj_(i, j); }

complex<double> Observables_MultiOrbSF::SiSjQ_Mean(int i, int j) { return SiSjQ_Mean_(i, j); }

complex<double> Observables_MultiOrbSF::SiSjQ_square_Mean(int i, int j) { return SiSjQ_square_Mean_(i, j); }

double Observables_MultiOrbSF::SiSj_Mean(int i, int j) { return SiSj_Mean_(i, j); }

double Observables_MultiOrbSF::SiSj_square_Mean(int i, int j) { return SiSj_square_Mean_(i, j); }

double Observables_MultiOrbSF::fermi_function(int n)
{
    double value;
    value = 1.0 / (exp(Parameters_.beta * (Hamiltonian_.eigs_[n] - Parameters_.mus)) + 1.0);
    return value;
}

void Observables_MultiOrbSF::calculate_quantum_SiSj()
{
    /*
    Matrix<complex<double>> F_u_u;
    Matrix<complex<double>> F_d_d;
    Matrix<complex<double>> F_u_d;
    Matrix<complex<double>> F_d_u;
    Matrix<complex<double>> omF_u_u;
    Matrix<complex<double>> omF_d_d;
    Matrix<complex<double>> omF_u_d;
    Matrix<complex<double>> omF_d_u;
    int nx, ny;
    int jx, jy;
    F_u_u.resize(ns_, ns_);
    F_d_d.resize(ns_, ns_);
    F_u_d.resize(ns_, ns_);
    F_d_u.resize(ns_, ns_);
    omF_u_u.resize(ns_, ns_);
    omF_d_d.resize(ns_, ns_);
    omF_u_d.resize(ns_, ns_);
    omF_d_u.resize(ns_, ns_);
    for (int i = 0; i < ns_; i++)
    {
        for (int j = 0; j < ns_; j++)
        {
            for (int n = 0; n < Hamiltonian_.eigs_.size(); n++)
            {
                F_u_u(i, j) += (conj(Hamiltonian_.Ham_(i, n)) * Hamiltonian_.Ham_(j, n) * fermi_function(n));
                F_d_d(i, j) += (conj(Hamiltonian_.Ham_(i + ns_, n)) * Hamiltonian_.Ham_(j + ns_, n) * fermi_function(n));
                F_u_d(i, j) += (conj(Hamiltonian_.Ham_(i, n)) * Hamiltonian_.Ham_(j + ns_, n) * fermi_function(n));
                F_d_u(i, j) += (conj(Hamiltonian_.Ham_(i + ns_, n)) * Hamiltonian_.Ham_(j, n) * fermi_function(n));
                omF_u_u(i, j) += (conj(Hamiltonian_.Ham_(i, n)) * Hamiltonian_.Ham_(j, n) * (1.0 - fermi_function(n)));
                omF_d_d(i, j) += (conj(Hamiltonian_.Ham_(i + ns_, n)) * Hamiltonian_.Ham_(j + ns_, n) * (1.0 - fermi_function(n)));
                omF_u_d(i, j) += (conj(Hamiltonian_.Ham_(i, n)) * Hamiltonian_.Ham_(j + ns_, n) * (1.0 - fermi_function(n)));
                omF_d_u(i, j) += (conj(Hamiltonian_.Ham_(i + ns_, n)) * Hamiltonian_.Ham_(j, n) * (1.0 - fermi_function(n)));
            }
        }
    }

    int i_;
    for (int ix = 0; ix < lx_; ix++)
    {
        for (int iy = 0; iy < ly_; iy++)
        {
            // i = Coordinates_.Nc(ix, iy);

            quantum_SiSj_(ix, iy) = 0.0;
            for (int j = 0; j < ns_; j++)
            {
                jx = Coordinates_.indx(j);
                jy = Coordinates_.indy(j);
                nx = (jx + ix) % lx_;
                ny = (jy + iy) % ly_;
                i_ = Coordinates_.Nc(nx, ny);
                quantum_SiSj_(ix, iy) += (
                0.25*(F_u_u(i_, i_) * F_u_u(j, j) + F_u_u(i_, j) * omF_u_u(j, i_)
                    - ( F_u_u(i_, i_) * F_d_d(j, j) + F_u_d(i_, j) * omF_d_u(j, i_) )
                    - ( F_d_d(i_, i_) * F_u_u(j, j) + F_d_u(i_, j) * omF_u_d(j, i_) )
                    + F_d_d(i_, i_) * F_d_d(j, j) + F_d_d(i_, j) * omF_d_d(j, i_))
                    + 0.5 * (F_u_d(i_, i_) * F_d_u(j, j) + F_u_u(i_, j) * omF_d_d(j, i_))
                    + 0.5 * (F_d_u(i_, i_) * F_u_d(j, j) + F_d_d(i_, j) * omF_u_u(j, i_))
                    ).real();
            }
            quantum_SiSj_(ix, iy) /= (ns_ * 1.0);
        }
    }

    //Fourier transform
    double phase, Cos_ij, Sin_ij;
    for (int qx = 0; qx < lx_; qx++)
    {
        for (int qy = 0; qy < ly_; qy++)
        {
            quantum_SiSjQ_(qx, qy) = complex<double>(0.0, 0.0);
            for (int xr = 0; xr < lx_; xr++)
            {
                for (int yr = 0; yr < ly_; yr++)
                {
                    phase = 2.0 * Parameters_.pi * (double(qx * xr) / double(lx_) + double(qy * yr) / double(ly_));
                    Cos_ij = cos(phase);
                    Sin_ij = sin(phase);
                    quantum_SiSjQ_(qx, qy) += quantum_SiSj_(xr, yr) * complex<double>(Cos_ij, Sin_ij);
                }
            }
            quantum_SiSjQ_(qx, qy) *= double(1.0 / (lx_ * ly_));
            //cout << qx << " "<< qy<< " "<<  SiSjQ_(qx,qy) << endl;
        }
    }

    */
}

void Observables_MultiOrbSF::calculate_local_density()
{
    int c1;
    complex<double> value = zero_complex;
    // cout <<"Parameter mus="<< Parameters_.mus<<endl;
    for (int i = 0; i < nbasis_; i++)
    {
        for (int sigma = 0; sigma < 2; sigma++)
        {
            local_density[i][sigma] = 0.0;
            c1 = i + (sigma * nbasis_);
            for (int n = 0; n < Hamiltonian_.eigs_.size(); n++)
            {
                local_density[i][sigma] += (conj(Hamiltonian_.Ham_(c1, n)) * Hamiltonian_.Ham_(c1, n) * fermi_function(n)).real();
            }

            // value += (conj(Hamiltonian_.Ham_(c1, 1)) * Hamiltonian_.Ham_(c1, 1));
        }
    }
}

void Observables_MultiOrbSF::SiSjFULL()
{


} // ----------

void Observables_MultiOrbSF::SiSjQ_Average()
{



} // ----------

void Observables_MultiOrbSF::quantum_SiSjQ_Average()
{


}

void Observables_MultiOrbSF::SiSj_Average()
{


} // ----------

void Observables_MultiOrbSF::quantum_SiSj_Average()
{


} // ----------

void Observables_MultiOrbSF::local_density_average()
{
    for (int i = 0; i < nbasis_; i++)
    {
        for (int sigma = 0; sigma < 2; sigma++)
        {
            local_density_Mean[i][sigma] += local_density[i][sigma];
            local_density_square_Mean[i][sigma] += pow(local_density[i][sigma], 2);
        }
    }
}

void Observables_MultiOrbSF::Total_Energy_Average(double Curr_QuantE, double CurrE)
{

    AVG_Total_Energy += Curr_QuantE + CurrE;
    AVG_Total_Energy_sqr += (Curr_QuantE + CurrE) * (Curr_QuantE + CurrE);
}

void Observables_MultiOrbSF::OccDensity(int tlabel)
{

} // ----------

void Observables_MultiOrbSF::DOSprint(int tlabel)
{

} // ----------

void Observables_MultiOrbSF::Initialize()
{

    complex<double> zero(0.0, 0.0);
    int space = 2 * ncells_ * n_orbs_;
    sx_.resize(space);
    sy_.resize(space);
    sz_.resize(space);

    local_density.resize(nbasis_);
    local_density_Mean.resize(nbasis_);
    local_density_square_Mean.resize(nbasis_);
    for (int i = 0; i < local_density.size(); i++)
    {
        local_density[i].resize(2);
        local_density_Mean[i].resize(2);
        local_density_square_Mean[i].resize(2);
    }

    Nematic_order_mean_ = 0.0;
    Nematic_order_square_mean_ = 0.0;

    SiSj_.resize(lx_, ly_);
    SiSj_Mean_.resize(lx_, ly_);
    SiSj_square_Mean_.resize(lx_, ly_);

    SiSjQ_Mean_.resize(lx_, ly_);
    SiSjQ_square_Mean_.resize(lx_, ly_);
    SiSjQ_.resize(lx_, ly_);

    quantum_SiSj_.resize(lx_, ly_);
    quantum_SiSj_Mean_.resize(lx_, ly_);
    quantum_SiSj_square_Mean_.resize(lx_, ly_);

    quantum_SiSjQ_Mean_.resize(lx_, ly_);
    quantum_SiSjQ_square_Mean_.resize(lx_, ly_);
    quantum_SiSjQ_.resize(lx_, ly_);

    for (int ix = 0; ix < lx_; ix++)
    {
        for (int iy = 0; iy < ly_; iy++)
        {
            SiSjQ_Mean_(ix, iy) = zero;
            SiSjQ_square_Mean_(ix, iy) = zero;
        }
    }

} // ----------

double Observables_MultiOrbSF::Omega(int i)
{
    return -20.0 + double(i) * dosincr_;
} // ----------


