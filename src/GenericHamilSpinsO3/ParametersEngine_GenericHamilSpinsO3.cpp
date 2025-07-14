#include "ParametersEngine_GenericHamilSpinsO3.h"

void Parameters_GenericHamilSpinsO3::Initialize(string inputfile_)
{

    cout << "____________________________________" << endl;
    cout << "Reading the inputfile: '" << inputfile_ <<"'"<< endl;
    cout << "____________________________________" << endl;


    ns = int(matchstring(inputfile_, "Nsites"));
    cout << "TotalNumberOfSites = " << ns << endl;


    lx = int(matchstring(inputfile_, "lx"));
    cout << "Lx = " << lx << endl;

    ly = int(matchstring(inputfile_, "ly"));
    cout << "Ly = " << ly << endl;


    PrintingNoOfTimeSlices=int(matchstring(inputfile_,"PrintingNoOfTimeSlices"));

    string genericconnectionsfiles_;
    genericconnectionsfiles_ = matchstring2(inputfile_, "GenericConnectionsFiles");
    stringstream genericconnectionsfiles_stream;
    genericconnectionsfiles_stream<<genericconnectionsfiles_;
    genericconnectionsfiles_stream>>N_ConnectionsFiles;
    cout<<"N_ConnectionsFiles = "<<N_ConnectionsFiles<<endl;
    ConnectionFiles.resize(N_ConnectionsFiles);
    string filename_temp;
    for(int i=0;i<N_ConnectionsFiles;i++){
        genericconnectionsfiles_stream>>filename_temp;
        ConnectionFiles[i]=filename_temp;
    }


    Connections.resize(N_ConnectionsFiles);
    for(int FileNo=0;FileNo<N_ConnectionsFiles;FileNo++){
        string line_connection;
        ifstream inputfileConnection(ConnectionFiles[FileNo].c_str());
        Connections[FileNo].clear();
        while(getline(inputfileConnection,line_connection)){
            Connections[FileNo].push_back(line_connection);
        }
    }


    cout << "____________________________________" << endl;
}

double Parameters_GenericHamilSpinsO3::matchstring(string file, string match)
{
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass = false;
    while (std::getline(readFile, line))
    {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass == false)
        {
            // ---------------------------------
            if (iss >> amount && test == match)
            {
                // cout << amount << endl;
                pass = true;
            }
            else
            {
                pass = false;
            }
            // ---------------------------------
            if (pass)
                break;
        }
    }
    if (pass == false)
    {
        string errorout = match;
        errorout += "= argument is missing in the input file!";
        throw std::invalid_argument(errorout);
    }
    cout << match << " = " << amount << endl;
    return amount;
}

string Parameters_GenericHamilSpinsO3::matchstring2(string file, string match)
{

    string line;
    ifstream readFile(file);
    string amount;
    int offset;

    if (readFile.is_open())
    {
        while (!readFile.eof())
        {
            getline(readFile, line);

            if ((offset = line.find(match, 0)) != string::npos)
            {
                amount = line.substr(offset + match.length() + 1);
            }
        }
        readFile.close();
    }
    else
    {
        cout << "Unable to open input file while in the Parameters class." << endl;
    }

    cout << match << " = " << amount << endl;
    return amount;
}


