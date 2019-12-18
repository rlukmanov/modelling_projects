#include "vmo.hpp"
#include "toposBuild.hpp"
#include <cstring>

using namespace std;

int main(int argc, char const *argv[])
{

    if (argc < 4) 
    {
        cout << "Solver v0.2 beta" << endl;
        cout << "Usage: <format> <topology> <threads>" << endl;
        cout << "Commands:" << endl;
        cout << "\t<format>                  Sparse matrix format. Choose between \"-csr\",\"-ellp\" and \"-coo\"." << endl; 
        cout << "\t<topology>                Topology which sparse matrix is generated from.Choose between \"-EE\" and \"-NN\"(EE is not supported yet)." << endl; 
        cout << "\t<threads>                 Threads to parallel SPMV operation." << endl;
        return 1;
    }
    else
    {
        if (atoi(argv[3]) < 1)
        {
            cout << "Threads number can't be lower than 1!" << endl;
            return 1;
        }

        size_t threads = atoi(argv[3]);
        if (!strcmp("-csr",argv[1]))
        {
            //csr
            cout << "Csr" <<endl;
            if (!strcmp("-EE",argv[2]))
            {
                //EE
                cout << "EE is not currently supported yet!" <<endl;
            }
            else  if (!strcmp("-NN",argv[2]))
            {
                //NN
                cout << "NN" <<endl;
            }
            else
            {
                cout << "Wrong topology chosen!" << endl;
                return 1;
            }
        }
        else if (!strcmp("-ellp",argv[1]))
        {
            //ellp
            cout << "Ellpack" <<endl;
            if (!strcmp("-EE",argv[2]))
            {
                //EE
                cout << "EE is not currently supported yet!" <<endl;
            }
            else  if (!strcmp("-NN",argv[2]))
            {
                //NN
                cout << "NN" <<endl;
            }
            else
            {
                cout << "Wrong topology chosen!" << endl;
                return 1;
            }
            
        }
        else if (!strcmp("-coo",argv[1]))
        {
            //coo
            cout << "Coord" <<endl;
            if (!strcmp("-EE",argv[2]))
            {
                //EE
                cout << "EE is not currently supported yet!" <<endl;
            }
            else  if (!strcmp("-NN",argv[2]))
            {
                //NN
                cout << "NN" <<endl;
            }
            else
            {
                cout << "Wrong topology chosen!" << endl;
                return 1;
            }
        }
        else
        {
            cout << "Wrong sparse matrix format chosen!" << endl;
            return 1;
        }

    }
}
