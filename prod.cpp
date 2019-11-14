#include <iostream>
#include <fstream>
#include <list>
#include <iterator>
#include <vector>
#include <string.h>
#include <omp.h>
#include <cstring>
#include <cstdlib>
#include "toposBuild.cpp"
#include "VariableSizeMeshContainer.cpp"
#include "FixedSizeMeshContainer.cpp"
#include "vtkGenerator.cpp"
#include "meshIO.cpp"

using namespace std;

int main(int argc, char **argv) {

    int Nx, Ny, k3, k4, nN, nE;
    double Lx, Ly;
    FixedSizeMeshContainer<double> C;
    VariableSizeMeshContainer<int> topoEN;
    VariableSizeMeshContainer<int> topoNE;
    VariableSizeMeshContainer<int> topoSN;
    VariableSizeMeshContainer<int> topoNS;
    int type = 0;

    double start, end;

    if ((argc == 1) || !((argc == 3) || (argc == 9) || (argc == 10))){
        cout << "Usage:\n\t-gen Lx Ly Nx Ny k3 k4 | -file\n\t-print (to stdout)| -out (to file) | -vtk \"filename\" " << endl;
        return 0;
    }
    try {
//creating mesh
        if (!strcmp(argv[1], "-gen")){
            Lx = atoi(argv[2]);
            Ly = atoi(argv[3]);
            Nx = atoi(argv[4]);
            Ny = atoi(argv[5]);
            k3 = atoi(argv[6]);
            k4 = atoi(argv[7]);
            if ((k3*k3 + k4*k4 == 0) || (k3 < 0) || (k4 < 0)) throw 1;

            C.setBlockSize(2);
            nN = Nx *  Ny;
            nE = num_elem(Nx, Ny, k3, k4);

            cout << "Time:" << endl;

            start = omp_get_wtime();
            topos::build_coord(C, Lx, Ly, Nx, Ny);
            end = omp_get_wtime();
            cout << "\tC: " << end - start << " sec" << endl;

            start = omp_get_wtime();
            topoEN = topos::build_topoEN(C, Nx, Ny, k3, k4, nE);
            end = omp_get_wtime();
            cout << "\ttopoEN: " << end - start << " sec" << endl;

            start = omp_get_wtime();
            topoNE = topos::build_topoNE(topoEN);
            end = omp_get_wtime();
            cout << "\ttopoNE: " << end - start << " sec" << endl;

            start = omp_get_wtime();
            topoSN = topos::build_topoSN(Nx,Ny);
            end = omp_get_wtime();
            cout << "\ttopoSN: " << end - start << " sec" << endl;

            start = omp_get_wtime();
            topoNS = topos::build_topoNS(topoSN);
            end = omp_get_wtime();
            cout << "\ttopoNS: " << end - start << " sec" << endl;
        }
        else if (!strcmp(argv[1], "-file")){
            if (read_file(C, topoEN,topoSN))
                throw 1;
            type = 1;
            topoNE = topos::build_topoNE(topoEN);
            topoNS = topos::build_topoNS(topoSN);
        }
        else {throw 1;}

        cout << "\nMemory:\n" << endl;
        cout << "\tC: " <<  sizeof(C) * sizeof(double) << " bytes" << endl;
        cout << "\ttopoEN: " <<  sizeof(topoEN) * sizeof(int) << " bytes" << endl;
        cout << "\ttopoNE: " <<  sizeof(topoNE) * sizeof(int) << " bytes" << endl;
        cout << "\ttopoSN: " <<  sizeof(topoSN) * sizeof(int) << " bytes" << endl;
        cout << "\ttopoNS: " <<  sizeof(topoNS) * sizeof(int) << " bytes" << endl;


//output mesh
        if (!strcmp(argv[argc-2], "-vtk"))
        {
            char* filename = argv[argc-1];
            strcat(filename,".vtk");
            vtkGenerator<double,int> vtk(filename);
            vtk.printInVTK(nN,C,topoEN);
        }
        else if (!strcmp(argv[argc-1], "-out")){
            write_file(C,topoEN,topoSN);
        }
        else if (!strcmp(argv[argc-1], "-print")){
            cout << "\nCoordinates:\n" << endl;
            C.printContainer();
            cout << "\nTopoEN:\n" << endl;
            topoEN.printContainer();
            cout << "\nTopoNE:\n" << endl;
            topoNE.printContainer();
            cout << "\nTopoSN:\n" << endl;
            topoSN.printContainer();
            cout << "\nTopoNS:\n" << endl;
            topoNS.printContainer();
            if (!type) draw_grid(Nx - 1, Ny - 1, k3, k4);
        }
        else {throw 1;}
    }
    catch(...){
        cout << "Error!" << endl;
        return 1;
    }
    
    return 0;
}
