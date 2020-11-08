#include "MFParams_MultiOrbSF.h"


void MFParams_MultiOrbSF::Adjust_MCWindow()
{
    double ratio;
    ratio = Parameters_.AccCount[0] / (Parameters_.AccCount[0] + Parameters_.AccCount[1]);
    //cout<<"ratio= "<< ratio << "temp= "<<Parameters_.temp << endl;
    Parameters_.AccCount[0] = 0;
    Parameters_.AccCount[1] = 0;
    Parameters_.WindowSize *= abs(1.0 + 1.0 * (ratio - 0.5));
    if(Parameters_.WindowSize>10){
        Parameters_.WindowSize=10.0;
    }
    //Parameters_.WindowSize =0.2;
    cout << "Ratio: " << ratio << "  window size:  " << Parameters_.WindowSize << endl;
    return;
} // ----------

void MFParams_MultiOrbSF::FieldThrow(int site, string mc_dof_type)
{
    int a, b;

    int Pi_multiple;

    double Pi = Parameters_.pi;
    double MC_Window = Parameters_.WindowSize;

    a = Coordinates_.indx_cellwise(site);
    b = Coordinates_.indy_cellwise(site);

    //ANGLES
    if (mc_dof_type == "phi")
    {
        ephi(a, b) += 2 * Pi * (random1() - 0.5) * MC_Window;

        Pi_multiple = ephi(a, b)/Pi;


        if (ephi(a, b) < 0.0)
        {
            ephi(a, b) = -ephi(a, b);
        }

        ephi(a, b) = fmod(ephi(a, b), 2.0 * Pi);


    }

    if (mc_dof_type == "theta")
    {
        etheta(a, b) += Pi * (random1() - 0.5) * MC_Window;
        if (etheta(a, b) < 0.0)
        {
            etheta(a, b) = -etheta(a, b);
        }

        etheta(a, b) = fmod(etheta(a, b),  Pi);

    }


    if (mc_dof_type == "theta_and_phi")
    {
        //phi
        ephi(a, b) += 2 * Pi * (random1() - 0.5) * MC_Window;

        Pi_multiple = ephi(a, b)/Pi;

        if (ephi(a, b) < 0.0)
        {
            ephi(a, b) = -ephi(a, b);
        }

        ephi(a, b) = fmod(ephi(a, b), 2.0 * Pi);


        //theta
        etheta(a, b) += Pi * (random1() - 0.5) * MC_Window;
        if (etheta(a, b) < 0.0)
        {
            etheta(a, b) = -etheta(a, b);
        }

        etheta(a, b) = fmod(etheta(a, b),  Pi);
    }


} // ----------

double MFParams_MultiOrbSF::random1()
{

    return dis1_(Generator1_);
}

double MFParams_MultiOrbSF::random2()
{

    return dis2_(Generator2_);
}

void MFParams_MultiOrbSF::initialize()
{

    lx_ = Coordinates_.lx_;
    ly_ = Coordinates_.ly_;

    // srand(Parameters_.RandomSeed);

    etheta_avg.resize(lx_, ly_);
    ephi_avg.resize(lx_, ly_);
    Disorder.resize(lx_, ly_);

    etheta.resize(lx_, ly_);
    ephi.resize(lx_, ly_);
    Moment_Size.resize(lx_, ly_);

    Sz.resize(lx_, ly_);
    Sx.resize(lx_, ly_);
    Sy.resize(lx_, ly_);

    ofstream Disorder_conf_file("Disorder_conf_used");
    Disorder_conf_file << "#seed=" << Parameters_.RandomDisorderSeed << " for mt19937_64 Generator is used" << endl;
    Disorder_conf_file << "#ix   iy    Dis[ix,iy]" << endl;

    ofstream Initial_MC_DOF_file("Initial_MC_DOF_values");

    Initial_MC_DOF_file << "#seed=" << Parameters_.RandomSeed << " for mt19937_64 Generator is used" << endl;
    Initial_MC_DOF_file << "#ix   iy   Theta(x,y)    Phi(x,y)      Moment_Size(x,y)" << endl;


    string temp_string;
    int ix_, iy_;
    if (Parameters_.Read_Seed_from_file_ == true)
    {
        ifstream Initial_Seed(Parameters_.Seed_file_name_);
        getline(Initial_Seed, temp_string);
        // cout << temp_string << endl;
        for (int ix = 0; ix < lx_; ix++)
        {
            for (int iy = 0; iy < ly_; iy++)
            {
                Initial_Seed >> ix_ >> iy_ >> etheta(ix, iy) >> ephi(ix, iy) >> Moment_Size(ix, iy);
                assert(ix_ == ix);
                assert(iy_ == iy);
            }
        }
    }

    else
    {
        for (int j = 0; j < ly_; j++)
        {
            for (int i = 0; i < lx_; i++)
            {
                //RANDOM fields
                if (Parameters_.MC_on_theta_and_phi == true)
                {
                    ephi(i, j) = 2.0 * random1() * PI;
                    etheta(i, j) = random1() * PI;
                }
                else
                {
                    if (Parameters_.MC_on_phi == true)
                    {
                        ephi(i, j) = 2.0 * random1() * PI;
                    }
                    else
                    {
                        ephi(i, j) = 0.0;
                    }

                    if (Parameters_.MC_on_theta == true)
                    {
                        etheta(i, j) = random1() * PI;
                    }
                    else
                    {
                        etheta(i, j) = 0.0;
                    }
                }


                Moment_Size(i, j) = 1.0;

                Initial_MC_DOF_file << i << setw(15) << j << setw(15) << etheta(i, j) << setw(15) << ephi(i, j)
                                    << setw(15) << Moment_Size(i, j) << endl;

            }
        }
    }

    //RANDOM Disorder
    for (int j = 0; j < ly_; j++)
    {
        for (int i = 0; i < lx_; i++)
        {
            Disorder(i, j) = Parameters_.Disorder_Strength * ((2.0 * random2()) - 1.0);
            Disorder_conf_file << i << "  " << j << "  " << Disorder(i, j) << endl;
        }
        Disorder_conf_file << endl;
    }

} // ----------

void MFParams_MultiOrbSF::Calculate_Fields_Avg()
{

    for (int j = 0; j < ly_; j++)
    {
        for (int i = 0; i < lx_; i++)
        {

            ephi_avg(i, j) = ephi_avg(i, j) + ephi(i, j);
            etheta_avg(i, j) = etheta_avg(i, j) + etheta(i, j);
        }
    }

} // ----------

void MFParams_MultiOrbSF::Read_classical_DOFs(string filename)
{

    string tmp_str;
    double tmp_double;
    ifstream fl_in(filename.c_str());
    getline (fl_in,tmp_str);

    for (int i = 0; i < lx_; i++)
    {
        for (int j = 0; j < ly_; j++)
        {
            fl_in >> tmp_double >> tmp_double >> etheta(i, j) >> ephi(i, j)>> tmp_double;
        }
    }

} // ----------

