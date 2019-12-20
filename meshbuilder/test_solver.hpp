#include <cstring>
#include <stdlib.h>
#include <omp.h>
#include "vmo.hpp"
#include "toposBuild.hpp"
#include "Sparse.hpp"
#include "IO.hpp"
#include <iostream>

using namespace std;

int solver(VariableSizeMeshContainer<int> topoNN, char* type_matr, int threads = 4) {
    double start, end;
    vector<double> x;
	ofstream fout;
    int n;

    if (threads < 1)
    {
        cout << "Threads number can't be lower than 1!" << endl;
        return 1;
    }

    if (!strcmp("csr", type_matr))
    {
        //csr
        SparseCSR<double> matrix(topoNN);

        n = matrix.getDenseRows();
        vector<double> b(n);

        matrix.set_model_values();

        cout << endl << "Type matrix: CSR\n" << endl;
        for (int i = 0; i < n; i++)
            b[i] = i+1;

        start = omp_get_wtime();
        x = vmo::conGradSolver(matrix,b,threads);
        end = omp_get_wtime();

        cout << "\nTime solution:\n\t" << end - start << " sec" << endl;

    }
    else if (!strcmp("ellp",type_matr))
    {
        //ellpack
        SparseELL<double> matrix(topoNN);

        n = matrix.getDenseRows();
        vector<double> b(n);

        matrix.set_model_values();

        cout << endl << "Type matrix: ELLPACK\n" << endl;
        for (int i = 0; i < n; i++)
            b[i] = i+1;

        start = omp_get_wtime();
        x = vmo::conGradSolver(matrix,b,threads);
        end = omp_get_wtime();

        cout << "\nTime solution:\n\t" << end - start << " sec" << endl;

    }
    else if (!strcmp("coo",type_matr))
    {
        //csr
        SparseCOO<double> matrix(topoNN);

        n = matrix.getDenseRows();
        vector<double> b(n);

        matrix.set_model_values();

        cout << endl << "Type matrix: COO\n" << endl;
        for (int i = 0; i < n; i++)
            b[i] = i+1;

        start = omp_get_wtime();
        x = vmo::conGradSolver(matrix,b,threads);
        end = omp_get_wtime();

        cout << "\nTime solution:\n\t" << end - start << " sec" << endl;

    }
    else
    {
        cout << "Wrong sparse matrix format chosen!" << endl;
        return 1;
    }

    fout.open("decision.txt");
    for (int i = 0; i < x.size(); i++)
		fout << x[i] << endl;
    fout.close();
	
	return 0;
}
