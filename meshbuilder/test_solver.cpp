#include <cstring>
#include <stdlib.h>
#include "vmo.hpp"
#include "toposBuild.hpp"
#include "Sparse.hpp"
#include "IO.hpp"

#define name_builder "prod"

using namespace std;

int main(int argc, char const *argv[])
{
    if (argc != 11)
    {
        cout << "Solver v0.2 beta" << endl;
        cout << "Usage: --set <Lx> <Ly> <Nx> <Ny> <k3> <k4> <format> <topology> <threads> | \n       --ex <Lx> <Ly> <Nx> <Ny> <k3> <k4> <format> <topology> <threads>"  << endl;
        cout << "Commands:" << endl;
        cout << "\t--ex                      An example of solving the  Ax = b using the conjugate gradient method. " << endl;
        cout << "\t--set                     The values of non-zero values and decision vectors are set by user" << endl;
        cout << "\t<format>                  Sparse matrix format. Choose between \"-csr\",\"-ellp\" and \"-coo\"." << endl;
        cout << "\t<topology>                Topology which sparse matrix is generated from.Choose between \"-EE\" and \"-NN\"." << endl;
        cout << "\t<threads>                 Threads to parallel SPMV operation." << endl;
        cout << "\t<Lx>                      Mesh X length" << endl;
        cout << "\t<Ly>                      Mesh Y length" << endl;
        cout << "\t<Nx>                      Number of X nodes" << endl;
        cout << "\t<Ny>                      Number of Y nodes" << endl;
        cout << "\t<k3>                      Number of squares in sequence" << endl;
        cout << "\t<k4>                      Number of triangles in sequence" << endl;
        return 1;
    }
    else
    {
        FixedSizeMeshContainer<double> C;//coordinates
        VariableSizeMeshContainer<int> topoEN;
        VariableSizeMeshContainer<int> topoSN;
        VariableSizeMeshContainer<int> topoNN;
        VariableSizeMeshContainer<int> topo;
        vector<double> x;
        int n;

        if (atoi(argv[10]) < 1)
        {
            cout << "Threads number can't be lower than 1!" << endl;
            return 1;
        }

        string command = "./";
        command += name_builder;
        command += " --gen ";
        for (int i = 2; i <= 7; i++) {
            command += argv[i];
            command += " ";
        }
        command += "--out";

        if (system(command.c_str()) != 0) return 1;
        cout << endl;

        if (read_file(C, topoEN,topoSN)) {
            cout << "File can't be read!" << endl;
            return 1;
        }

        size_t threads = atoi(argv[10]);
        if (!strcmp("-csr",argv[8]))
        {
            //csr
            if (!strcmp("-EE",argv[9]))
            {
                //EE
                cout << "EE is not currently supported yet!" <<endl;
                return 1;
            }
            else  if (!strcmp("-NN",argv[9]))
            {
                //NN
                topoNN = topos::build_topoNN(topoSN);
                SparseCSR<double> matrix(topoNN);

                if (!strcmp(argv[1],"--set")){
                    double elem;
                    vector<double> a(matrix.getValuesSize());
                    vector<double> b;
                    n = matrix.getDenseRows();

                    M1:
                    cout << "Enter a vector of nonzero matrix elements (" << matrix.getValuesSize() << " elements):" << endl;
                    for (size_t i = 0; i < matrix.getValuesSize(); i++){
                        cout << "\ta[" << i << "] = ";
                        cin >> a[i];
                    }
                    cout << endl;
                    matrix.setValues(a);
                    cout << "CSR" << endl;
                    matrix.outputMatrix();
                    if (!vmo::isPositive(matrix.getDenseMatrix(),n)){
                        cout << "Unfortunately, the Sylvester criterion is not fulfilled. Re-enter nonzero matrix elements" << endl << endl;
                        goto M1;
                    }

                    cout << "Enter the resulting vector:" << endl;
                    for (int i = 0; i < n; i++){
                        cout << "\tb[" << i << "] = ";
                        cin >> elem;
                        b.push_back(elem);
                    }
                    cout << endl;

                    x = vmo::conGradSolver(matrix,b,threads);
                }
                else if (!strcmp(argv[1],"--ex")){
                    n = matrix.getDenseRows();
                    vector<double> a(matrix.getValuesSize());
                    vector<double> b(n);

                    for (size_t i = 0; i < matrix.getValuesSize(); i++)
                        a[i] = 1;
                    a[0] = a[matrix.getValuesSize()-1] = 5;
                    matrix.setValues(a);
                    cout << "CSR" << endl;
                    matrix.outputMatrix();
                    for (int i = 0; i < n; i++)
                        b[i] = i+1;

                    x = vmo::conGradSolver(matrix,b,threads);
                }
                else {
                    cout << "Unknown first argument!" << endl;
                    return 1;
                }
            }
            else
            {
                cout << "Wrong topology chosen!" << endl;
                return 1;
            }
        }
        else if (!strcmp("-ellp",argv[8]))
        {
            //ellp
            if (!strcmp("-EE",argv[9]))
            {
                //EE
                cout << "EE is not currently supported yet!" <<endl;
                return 1;
            }
            else  if (!strcmp("-NN",argv[9]))
            {
                //NN
                topoNN = topos::build_topoNN(topoSN);
                SparseELL<double> matrix(topoNN);

                if (!strcmp(argv[1],"--set")){
                    double elem;
                    vector<double> a(matrix.getValuesSize());
                    vector<double> b;
                    n = matrix.getDenseRows();

                    M2:
                    cout << "Enter a vector of nonzero matrix elements (" << matrix.getValuesSize() << " elements):" << endl;
                    for (size_t i = 0; i < matrix.getValuesSize(); i++){
                        cout << "\ta[" << i << "] = ";
                        cin >> a[i];
                    }
                    cout << endl;
                    matrix.setValues(a);
                    cout << "Ellpack" <<endl;
                    matrix.outputMatrix();
                    if (!vmo::isPositive(matrix.getDenseMatrix(),n)){
                        cout << "Unfortunately, the Sylvester criterion is not fulfilled. Re-enter nonzero matrix elements" << endl << endl;
                        goto M2;
                    }

                    cout << "Enter the resulting vector:" << endl;
                    for (int i = 0; i < n; i++){
                        cout << "\tb[" << i << "] = ";
                        cin >> elem;
                        b.push_back(elem);
                    }
                    cout << endl;

                    x = vmo::conGradSolver(matrix,b,threads);
                }
                else if (!strcmp(argv[1],"--ex")){
                    n = matrix.getDenseRows();
                    vector<double> a(matrix.getValuesSize());
                    vector<double> b(n);

                    for (size_t i = 0; i < matrix.getValuesSize(); i++)
                        a[i] = 1;
                    a[0] = a[matrix.getValuesSize()-1] = 5;
                    matrix.setValues(a);
                    cout << "Ellpack" <<endl;
                    matrix.outputMatrix();
                    for (int i = 0; i < n; i++)
                        b[i] = i+1;

                    x = vmo::conGradSolver(matrix,b,threads);
                }
                else {
                    cout << "Unknown first argument!" << endl;
                    return 1;
                }
            }
            else
            {
                cout << "Wrong topology chosen!" << endl;
                return 1;
            }
            
        }
        else if (!strcmp("-coo",argv[8]))
        {
            //coo
            if (!strcmp("-EE",argv[9]))
            {
                //EE
                cout << "EE is not currently supported yet!" <<endl;
                return 1;
            }
            else  if (!strcmp("-NN",argv[9]))
            {
                //NN
                topoNN = topos::build_topoNN(topoSN);
                SparseCOO<double> matrix(topoNN);

                if (!strcmp(argv[1],"--set")){
                    double elem;
                    vector<double> a(matrix.getValuesSize());
                    vector<double> b;
                    n = matrix.getDenseRows();

                    M3:
                    cout << "Enter a vector of nonzero matrix elements (" << matrix.getValuesSize() << " elements):" << endl;
                    for (size_t i = 0; i < matrix.getValuesSize(); i++){
                        cout << "\ta[" << i << "] = ";
                        cin >> a[i];
                    }
                    cout << endl;
                    matrix.setValues(a);
                    cout << "Coord" <<endl;
                    matrix.outputMatrix();
                    if (!vmo::isPositive(matrix.getDenseMatrix(),n)){
                        cout << "Unfortunately, the Sylvester criterion is not fulfilled. Re-enter nonzero matrix elements" << endl << endl;
                        goto M3;
                    }

                    cout << "Enter the resulting vector:" << endl;
                    for (int i = 0; i < n; i++){
                        cout << "\tb[" << i << "] = ";
                        cin >> elem;
                        b.push_back(elem);
                    }
                    cout << endl;

                    x = vmo::conGradSolver(matrix,b,threads);
                }
                else if (!strcmp(argv[1],"--ex")){
                    n = matrix.getDenseRows();
                    vector<double> a(matrix.getValuesSize());
                    vector<double> b(n);

                    for (size_t i = 0; i < matrix.getValuesSize(); i++)
                        a[i] = 1;
                    a[0] = a[matrix.getValuesSize()-1] = 5;
                    matrix.setValues(a);
                    cout << "Coord" <<endl;
                    matrix.outputMatrix();
                    for (int i = 0; i < n; i++)
                        b[i] = i+1;

                    x = vmo::conGradSolver(matrix,b,threads);
                }
                else {
                    cout << "Unknown first argument!" << endl;
                    return 1;
                }
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


        cout << endl << "Solutions vector:" << endl;
        for (int i = 0; i < n; i++)
            cout << "\tx["<< i << "] = " << x[i] << endl;
        cout << endl;
    }
}
