#include "Hamiltonian_MultiOrbSF.h"

double Hamiltonian_MultiOrbSF::chemicalpotential(double muin, double filling)
{
    double mu_out;
    double n1, N;
    double dMubydN;
    double nstate = eigs_.size();
    dMubydN = 0.05 * (eigs_[nstate - 1] - eigs_[0]) / nstate;
    N = filling * double(eigs_.size());
    //temp=Parameters_.temp;
    mu_out = muin;
    bool converged = false;

    if (!Parameters_.fix_mu)
    {
        assert(!Parameters_.fix_mu);

        if (1 == 2)
        {
            for (int i = 0; i < 100000; i++)
            {
                n1 = 0.0;
                for (int j = 0; j < nstate; j++)
                {
                    n1 += double(1.0 / (exp((eigs_[j] - mu_out) * Parameters_.beta) + 1.0));
                }
                //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
                if (abs(N - n1) < double(0.00001))
                {
                    //cout<<abs(N-n1)<<endl;
                    converged = true;
                    break;
                }
                else
                {
                    mu_out += (N - n1) * dMubydN;
                    //cout<<i<<"    "<<n1<<"    "<<N-n1<<endl;
                }
            }

            if (!converged)
            {
                cout << "mu_not_converged, N = " << n1 << endl;
            }
            else
            {
                //cout<<"mu converged, N = "<<n1<<endl;
            }
        }

        double mu1, mu2;
        double mu_temp = muin;
        //cout<<"mu_input = "<<mu_temp<<endl;
        if (1 == 1)
        {
            mu1 = eigs_[0];
            mu2 = eigs_[nstate - 1];
            for (int i = 0; i < 4000000; i++)
            {
                n1 = 0.0;
                for (int j = 0; j < nstate; j++)
                {
                    n1 += double(1.0 / (exp((eigs_[j] - mu_temp) * Parameters_.beta) + 1.0));
                }
                //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
                if (abs(N - n1) < double(0.00001))
                {
                    //cout<<abs(N-n1)<<endl;
                    converged = true;
                    break;
                }
                else
                {
                    if (n1 > N)
                    {
                        mu2 = mu_temp;
                        mu_temp = 0.5 * (mu1 + mu_temp);
                    }
                    else
                    {
                        mu1 = mu_temp;
                        mu_temp = 0.5 * (mu2 + mu_temp);
                    }
                }
                //cout<<"mu_temp = "<<mu_temp<<"   "<<n1<<endl;
            }

            if (!converged)
            {
                cout << "mu_not_converged, N = " << n1 << endl;
            }
            else
            {
                //cout<<"mu converged, N = "<<n1<<endl;
            }

            mu_out = mu_temp;
        }

        return mu_out;
    }
    else
    {
        assert(Parameters_.fix_mu);
        return Parameters_.fixed_mu_value;
    }
} // ----------

double Hamiltonian_MultiOrbSF::chemicalpotentialCluster(double muin, double filling)
{
    double mu_out;
    double n1, N;
    double dMubydN;
    double nstate = eigsCluster_.size();
    dMubydN = 0.05 * (eigsCluster_[nstate - 1] - eigsCluster_[0]) / nstate;
    N = filling * double(eigsCluster_.size());
    //temp=Parameters_.temp;
    mu_out = muin;
    bool converged = false;

    if (!Parameters_.fix_mu)
    {
        assert(!Parameters_.fix_mu);
        if (1 == 2)
        {
            for (int i = 0; i < 100000; i++)
            {
                n1 = 0.0;
                for (int j = 0; j < nstate; j++)
                {
                    n1 += double(1.0 / (exp((eigsCluster_[j] - mu_out) * Parameters_.beta) + 1.0));
                }
                //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
                if (abs(N - n1) < double(0.00001))
                {
                    //cout<<abs(N-n1)<<endl;
                    converged = true;
                    break;
                }
                else
                {
                    mu_out += (N - n1) * dMubydN;
                    //cout<<i<<"    "<<n1<<"    "<<N-n1<<endl;
                }
            }

            if (!converged)
            {
                cout << "mu_not_converged, N = " << n1 << endl;
            }
            else
            {
                //cout<<"mu converged, N = "<<n1<<endl;
            }
        }

        double mu1, mu2;
        double mu_temp = muin;
        //cout<<"mu_input = "<<mu_temp<<endl;
        if (1 == 1)
        {
            mu1 = eigsCluster_[0];
            mu2 = eigsCluster_[nstate - 1];
            for (int i = 0; i < 4000000; i++)
            {
                n1 = 0.0;
                for (int j = 0; j < nstate; j++)
                {
                    n1 += double(1.0 / (exp((eigsCluster_[j] - mu_temp) * Parameters_.beta) + 1.0));
                }
                //cout <<"i  "<< i << "  n1  " << n1 << "  mu  " << mu_out<< endl;
                if (abs(N - n1) < double(0.00001))
                {
                    //cout<<abs(N-n1)<<endl;
                    converged = true;
                    break;
                }
                else
                {
                    if (n1 > N)
                    {
                        mu2 = mu_temp;
                        mu_temp = 0.5 * (mu1 + mu_temp);
                    }
                    else
                    {
                        mu1 = mu_temp;
                        mu_temp = 0.5 * (mu2 + mu_temp);
                    }
                }
                //cout<<"mu_temp = "<<mu_temp<<"   "<<n1<<endl;
            }

            if (!converged)
            {
                cout << "mu_not_converged, N = " << n1 << endl;
            }
            else
            {
                //cout<<"mu converged, N = "<<n1<<endl;
            }

            mu_out = mu_temp;
        }

        return mu_out;
    }
    else
    {
        assert(Parameters_.fix_mu);
        return Parameters_.fixed_mu_value;
    }
} // ----------

void Hamiltonian_MultiOrbSF::Initialize()
{

    ly_cluster_ = Parameters_.ly_cluster;
    lx_cluster_ = Parameters_.lx_cluster;
    ncells_cluster = lx_cluster_*ly_cluster_;

    ly_ = Parameters_.ly;
    lx_ = Parameters_.lx;
    ncells_ = lx_*ly_;
    n_orbs_ = Parameters_.n_orbs;
    int space = 2 * ncells_ * n_orbs_;
    int spaceCluster = 2 * ncells_cluster* n_orbs_;

    HTB_.resize(space, space);
    Ham_.resize(space, space);
    HTBCluster_.resize(spaceCluster, spaceCluster);
    HamCluster_.resize(spaceCluster, spaceCluster);
    eigs_.resize(space);
    sx_.resize(space);
    sy_.resize(space);
    sz_.resize(space);
    eigs_saved_.resize(space);
    eigsCluster_.resize(spaceCluster);
    eigsCluster_saved_.resize(spaceCluster);

} // ----------

double Hamiltonian_MultiOrbSF::TotalDensity()
{
    double n1 = 0.0;
    for (int j = 0; j < eigs_.size(); j++)
    {
        n1 += 1.0 / (exp(Parameters_.beta * (eigs_[j] - Parameters_.mus * 1.0)) + 1.0);
    }
    return n1;
} // ----------

double Hamiltonian_MultiOrbSF::ClusterDensity()
{
    double n1 = 0.0;
    for (int j = 0; j < eigsCluster_.size(); j++)
    {
        n1 += 1.0 / (exp(Parameters_.beta * (eigsCluster_[j] - Parameters_.mus_Cluster * 1.0)) + 1.0);
    }
    return n1;
} // ----------

double Hamiltonian_MultiOrbSF::E_QM()
{
    double E = 0.0;
    for (int j = 0; j < eigs_.size(); j++)
    {
        E += (eigs_[j]) / (exp(Parameters_.beta * (eigs_[j] - Parameters_.mus)) + 1.0);
    }
    return E;
} // ----------

double Hamiltonian_MultiOrbSF::E_QMCluster()
{
    double E = 0.0;
    for (int j = 0; j < eigsCluster_.size(); j++)
    {
        E += (eigsCluster_[j]) / (exp(Parameters_.beta * (eigsCluster_[j] - Parameters_.mus_Cluster)) + 1.0);
    }
    return E;
} // ----------

double Hamiltonian_MultiOrbSF::GetCLEnergy()
{

    double EClassical;
    int cell;
    double ei, ai;

    for (int i = 0; i < lx_; i++)
    {
        for (int j = 0; j < ly_; j++)
        {
            cell = Coordinates_.Ncell(i, j);
            ei = MFParams_.etheta(i, j);
            ai = MFParams_.ephi(i, j);
            sx_[cell] = MFParams_.Moment_Size(i, j) * cos(ai) * sin(ei);
            sy_[cell] = MFParams_.Moment_Size(i, j) * sin(ai) * sin(ei);
            sz_[cell] = MFParams_.Moment_Size(i, j) * cos(ei);
        }
    }

    // Classical Energy
    EClassical = double(0.0);

    int _ix, _iy;
    for (int i = 0; i < ncells_; i++)
    {
        _ix = Coordinates_.indx_cellwise(i);
        _iy = Coordinates_.indy_cellwise(i);

        cell = Coordinates_.neigh(i, 0); //+x
        EClassical += 1.0 * Parameters_.K1x * ( (sx_[i] * sx_[cell]) + (sy_[i] * sy_[cell]) + (1.0 * sz_[i] * sz_[cell]));
        cell = Coordinates_.neigh(i, 2); //+y
        EClassical += Parameters_.K1y * ((sx_[i] * sx_[cell]) + (sy_[i] * sy_[cell]) + (1.0 * sz_[i] * sz_[cell]));
    }

    return EClassical;
} // ----------

void Hamiltonian_MultiOrbSF::InteractionsCreate()
{

    int a;
    double ei, ai;
    int index;
    int i_posx, i_posy;

    Ham_ = HTB_;
    // Ham_.print();

    for (int i = 0; i < ncells_; i++)
    { // For each cell
        i_posx = Coordinates_.indx_cellwise(i);
        i_posy = Coordinates_.indy_cellwise(i);
        ei = MFParams_.etheta(i_posx, i_posy);
        ai = MFParams_.ephi(i_posx, i_posy);

        for(int orb=0;orb<n_orbs_;orb++){

            index=Coordinates_.Nbasis(i_posx, i_posy, orb);

            Ham_(index, index) += Parameters_.J_Hund[orb] * (cos(ei)) * 0.5 * MFParams_.Moment_Size(i_posx, i_posy);
            Ham_(index + (ncells_*n_orbs_), index + (ncells_*n_orbs_)) += Parameters_.J_Hund[orb] * (-cos(ei)) * 0.5 * MFParams_.Moment_Size(i_posx, i_posy);
            Ham_(index, index + (ncells_*n_orbs_)) += Parameters_.J_Hund[orb] * sin(ei) * complex<double>(cos(ai), -sin(ai)) * 0.5 * MFParams_.Moment_Size(i_posx, i_posy); //S-
            Ham_(index + (ncells_*n_orbs_), index) += Parameters_.J_Hund[orb] * sin(ei) * complex<double>(cos(ai), sin(ai)) * 0.5 * MFParams_.Moment_Size(i_posx, i_posy);  //S+

            // On-Site potential
            for (int spin = 0; spin < 2; spin++)
            {
                a = Coordinates_.Nbasis(i_posx,i_posy,orb) + ncells_*n_orbs_*spin;
                Ham_(a, a) += complex<double>(1.0, 0.0) * (
                            Parameters_.OnSiteE[orb] +
                            MFParams_.Disorder(i_posx, i_posy)
                            );
            }
        }
    }

} // ----------

void Hamiltonian_MultiOrbSF::InteractionsClusterCreate(int Center_site)
{

    int x_pos, y_pos;
    double ei, ai;
    int a;
    int i_original;
    int index;
    int i_posx, i_posy;

    HamCluster_ = HTBCluster_;

    for (int i = 0; i < ncells_cluster; i++)
    { // For each cell in cluster

        i_posx = CoordinatesCluster_.indx_cellwise(i);
        i_posy = CoordinatesCluster_.indy_cellwise(i);

        x_pos = Coordinates_.indx_cellwise(Center_site) - int(Parameters_.lx_cluster / 2) + CoordinatesCluster_.indx_cellwise(i);
        y_pos = Coordinates_.indy_cellwise(Center_site) - int(Parameters_.ly_cluster / 2) + CoordinatesCluster_.indy_cellwise(i);
        x_pos = (x_pos + Coordinates_.lx_) % Coordinates_.lx_;
        y_pos = (y_pos + Coordinates_.ly_) % Coordinates_.ly_;

        i_original=Coordinates_.Ncell(x_pos, y_pos);
        ei = MFParams_.etheta(x_pos, y_pos);
        ai = MFParams_.ephi(x_pos, y_pos);

        for(int orb=0;orb<n_orbs_;orb++){

            index=CoordinatesCluster_.Nbasis(i_posx, i_posy, orb);

            HamCluster_(index, index) += Parameters_.J_Hund[orb] * (cos(ei)) * 0.5 * MFParams_.Moment_Size(x_pos, y_pos);
            HamCluster_(index + ncells_cluster*n_orbs_, index + ncells_cluster*n_orbs_) += Parameters_.J_Hund[orb] * (-cos(ei)) * 0.5 * MFParams_.Moment_Size(x_pos, y_pos);
            HamCluster_(index, index + ncells_cluster*n_orbs_) += Parameters_.J_Hund[orb] * sin(ei) * complex<double>(cos(ai), -sin(ai)) * 0.5 * MFParams_.Moment_Size(x_pos, y_pos); //S-
            HamCluster_(index + ncells_cluster*n_orbs_, index) += Parameters_.J_Hund[orb] * sin(ei) * complex<double>(cos(ai), sin(ai)) * 0.5 * MFParams_.Moment_Size(x_pos, y_pos);  //S+


            for (int spin = 0; spin < 2; spin++)
            {
                a = CoordinatesCluster_.Nbasis(i_posx,i_posy,orb) + ncells_cluster*n_orbs_*spin;
                HamCluster_(a, a) += complex<double>(1.0, 0.0) * (
                            Parameters_.OnSiteE[orb] +
                            MFParams_.Disorder(x_pos, y_pos)
                            );
            }
        }
    }


} // ----------

void Hamiltonian_MultiOrbSF::Check_up_down_symmetry()

{
    complex<double> temp(0, 0);
    complex<double> temp2;

    for (int i = 0; i < ncells_*n_orbs_; i++)
    {
        for (int j = 0; j < ncells_*n_orbs_; j++)
        {
            temp2 = Ham_(i, j) - Ham_(i + ncells_*n_orbs_, j + ncells_*n_orbs_); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            temp += temp2 * conj(temp2);
        }
    }

    cout << "Assymetry in up-down sector: " << temp << endl;
}

void Hamiltonian_MultiOrbSF::Check_Hermiticity()

{
    complex<double> temp(0, 0);
    complex<double> temp2;

    for (int i = 0; i < HamCluster_.n_row(); i++)
    {
        for (int j = 0; j < HamCluster_.n_row(); j++)
        {
            if (HamCluster_(i, j) != conj(HamCluster_(j, i)))
            {
                cout << i << "," << j << endl;
                cout << "i,j = " << HamCluster_(i, j) << ", j,i=" << conj(HamCluster_(j, i)) << endl;
            }
            assert(HamCluster_(i, j) == conj(HamCluster_(j, i))); //+ Ham_(i+orbs_*ns_,j) + Ham_(i,j+orbs_*ns_);
            //temp +=temp2*conj(temp2);
        }
    }

    // cout<<"Hermiticity: "<<temp<<endl;
}

void Hamiltonian_MultiOrbSF::Diagonalize(char option)
{

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);

    char jobz = option;
    // jobz = 'V';
    char uplo = 'U'; //WHY ONLY 'L' WORKS?
    int n = Ham_.n_row();
    int lda = Ham_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3 * n - 2);
    int info;
    int lwork = -1;

    eigs_.resize(Ham_.n_row());
    fill(eigs_.begin(), eigs_.end(), 0);
    // query:
    zheev_(&jobz, &uplo, &n, &(Ham_(0, 0)), &lda, &(eigs_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz, &uplo, &n, &(Ham_(0, 0)), &lda, &(eigs_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    if (info != 0)
    {
        std::cerr << "info=" << info << "\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}
}

void Hamiltonian_MultiOrbSF::DiagonalizeCluster(char option)
{

    //extern "C" void   zheev_(char *,char *,int *,std::complex<double> *, int *, double *,
    //                       std::complex<double> *,int *, double *, int *);

    char jobz = option;
    // jobz = 'V';
    // cout << option;
    char uplo = 'U'; //WHY ONLY 'L' WORKS?
    int n = HamCluster_.n_row();
    int lda = HamCluster_.n_col();
    vector<complex<double>> work(3);
    vector<double> rwork(3 * n - 2);
    int info;
    int lwork = -1;

    eigsCluster_.resize(HamCluster_.n_row());
    fill(eigsCluster_.begin(), eigsCluster_.end(), 0);
    // query:
    zheev_(&jobz, &uplo, &n, &(HamCluster_(0, 0)), &lda, &(eigsCluster_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    //lwork = int(real(work[0]))+1;
    lwork = int((work[0].real()));
    work.resize(lwork);
    // real work:
    zheev_(&jobz, &uplo, &n, &(HamCluster_(0, 0)), &lda, &(eigsCluster_[0]), &(work[0]), &lwork, &(rwork[0]), &info);
    if (info != 0)
    {
        std::cerr << "info=" << info << "\n";
        perror("diag: zheev: failed with info!=0.\n");
    }

    // Ham_.print();

    //  for(int i=0;i<eigs_.size();i++){
    //    cout<<eigs_[i]<<endl;
    //}
}

void Hamiltonian_MultiOrbSF::HTBCreate()
{

    //Convention used
    //orb=0=d
    //orb=1=px
    //orb=2=py

    int mx = Parameters_.TBC_mx;
    int my = Parameters_.TBC_my;

    complex<double> phasex, phasey;
    int l, m, a, b;
    int lx_pos, ly_pos;
    int mx_pos, my_pos;

    HTB_.fill(0.0);

    for (l = 0; l < ncells_; l++)
    {
        lx_pos = Coordinates_.indx_cellwise(l);
        ly_pos = Coordinates_.indy_cellwise(l);


        // * +x direction Neighbor
        if (lx_pos == (Coordinates_.lx_ - 1))
        {
            phasex = Parameters_.BoundaryConnection*exp(iota_complex * 2.0 * (1.0 * mx) * PI / (1.0 * Parameters_.TBC_cellsX));
            phasey = one_complex;
        }
        else
        {
            phasex = one_complex;
            phasey = one_complex;
        }
        m = Coordinates_.neigh(l, 0); //+x neighbour cell
        mx_pos = Coordinates_.indx_cellwise(m);
        my_pos = Coordinates_.indy_cellwise(m);

        for (int spin=0; spin<2; spin++){
            for(int orb1=0;orb1<n_orbs_;orb1++){
                for(int orb2=0;orb2<n_orbs_;orb2++){
                    if(Parameters_.hopping_NN_X(orb1,orb2)!=0.0){
                        a = Coordinates_.Nbasis(lx_pos,ly_pos,orb1) + ncells_*n_orbs_*spin;
                        b = Coordinates_.Nbasis(mx_pos,my_pos,orb2) + ncells_*n_orbs_*spin;
                        assert(a != b);
                        if (a != b)
                        {
                            HTB_(b, a) = complex<double>(1.0 *Parameters_.hopping_NN_X(orb1,orb2), 0.0) * phasex;
                            HTB_(a, b) = conj(HTB_(b, a));
                        }
                    }
                }
            }
        }

        // * +y direction Neighbor
        if (ly_pos == (Coordinates_.ly_ - 1))
        {
            phasex = one_complex;
            phasey = Parameters_.BoundaryConnection*exp(iota_complex * 2.0 * (1.0 * my) * PI / (1.0 * Parameters_.TBC_cellsY));
        }
        else
        {
            phasex = one_complex;
            phasey = one_complex;
        }
        m = Coordinates_.neigh(l, 2); //+y neighbour cell
        mx_pos = Coordinates_.indx_cellwise(m);
        my_pos = Coordinates_.indy_cellwise(m);

        for (int spin = 0; spin < 2; spin++){
            for(int orb1=0;orb1<n_orbs_;orb1++){
                for(int orb2=0;orb2<n_orbs_;orb2++){
                    if(Parameters_.hopping_NN_Y(orb1,orb2)!=0.0){

                        a = Coordinates_.Nbasis(lx_pos,ly_pos,orb1) + ncells_*n_orbs_*spin;
                        b = Coordinates_.Nbasis(mx_pos,my_pos,orb2) + ncells_*n_orbs_*spin;
                        assert(a != b);
                        if (a != b)
                        {
                            HTB_(b, a) = complex<double>(1.0*Parameters_.hopping_NN_Y(orb1,orb2), 0.0) * phasey;
                            HTB_(a, b) = conj(HTB_(b, a));
                        }
                    }
                }
            }
        }



    }




} // ----------

void Hamiltonian_MultiOrbSF::HTBClusterCreate()
{

    if(Parameters_.ED_==false){

        int l, m, a, b;
        int lx_pos, ly_pos;
        int mx_pos, my_pos;

        HTBCluster_.fill(0.0);

        for (l = 0; l < ncells_cluster; l++)
        {

            lx_pos = CoordinatesCluster_.indx_cellwise(l);
            ly_pos = CoordinatesCluster_.indy_cellwise(l);


            // * +x direction Neighbor
            m = CoordinatesCluster_.neigh(l, 0);
            mx_pos = CoordinatesCluster_.indx_cellwise(m);
            my_pos = CoordinatesCluster_.indy_cellwise(m);

            for (int spin = 0; spin < 2; spin++){
                for(int orb1=0;orb1<n_orbs_;orb1++){
                    for(int orb2=0;orb2<n_orbs_;orb2++){
                        if(Parameters_.hopping_NN_X(orb1,orb2)!=0.0){
                            a = CoordinatesCluster_.Nbasis(lx_pos,ly_pos,orb1) + ncells_cluster*n_orbs_*spin;
                            b = CoordinatesCluster_.Nbasis(mx_pos,my_pos,orb2) + ncells_cluster*n_orbs_*spin;

                            assert(a != b);
                            if (a != b)
                            {
                                HTBCluster_(b, a) = complex<double>(1.0 *Parameters_.hopping_NN_X(orb1,orb2), 0.0);
                                HTBCluster_(a, b) = conj(HTBCluster_(b, a));
                            }
                        }
                    }
                }
            }

            // * +y direction Neighbor
            m = CoordinatesCluster_.neigh(l, 2);
            mx_pos = CoordinatesCluster_.indx_cellwise(m);
            my_pos = CoordinatesCluster_.indy_cellwise(m);

            for (int spin = 0; spin < 2; spin++){
                for(int orb1=0;orb1<n_orbs_;orb1++){
                    for(int orb2=0;orb2<n_orbs_;orb2++){
                        if(Parameters_.hopping_NN_Y(orb1,orb2)!=0.0){
                            a = CoordinatesCluster_.Nbasis(lx_pos,ly_pos,orb1) + ncells_cluster*n_orbs_*spin;
                            b = CoordinatesCluster_.Nbasis(mx_pos,my_pos,orb2) + ncells_cluster*n_orbs_*spin;
                            assert(a != b);
                            if (a != b)
                            {
                                HTBCluster_(b, a) = complex<double>(1.0 *Parameters_.hopping_NN_Y(orb1,orb2), 0.0);
                                HTBCluster_(a, b) = conj(HTBCluster_(b, a));
                            }
                        }
                    }
                }
            }

        }

    }
    else{
        HTBCluster_=HTB_;
    }

    // HTBCluster_.print();

} // ----------

void Hamiltonian_MultiOrbSF::Hoppings()
{
//Using matrices from Parameters_

} // ----------

void Hamiltonian_MultiOrbSF::copy_eigs(int i)
{

    int space = 2 * ncells_ *n_orbs_;

    if (i == 0)
    {
        for (int j = 0; j < space; j++)
        {
            eigs_[j] = eigs_saved_[j];
        }
    }
    else
    {
        for (int j = 0; j < space; j++)
        {
            eigs_saved_[j] = eigs_[j];
        }
    }
}

void Hamiltonian_MultiOrbSF::copy_eigs_Cluster(int i)
{

    int space = 2 * ncells_cluster *n_orbs_;

    if (i == 0)
    {
        for (int j = 0; j < space; j++)
        {
            eigsCluster_[j] = eigsCluster_saved_[j];
        }
    }
    else
    {
        for (int j = 0; j < space; j++)
        {
            eigsCluster_saved_[j] = eigsCluster_[j];
        }
    }
}


