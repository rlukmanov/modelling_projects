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
#include "IO.cpp"

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

    if ((argc == 1) || !((argc == 3) || (argc == 9) || (argc == 10)) || (strcmp(argv[argc-1], "--print") && strcmp(argv[argc-1], "--out") && strcmp(argv[argc-2], "--out") && strcmp(argv[argc-2], "--vtk"))){
        cout << "Mesh builder v0.1 beta" << endl;
        cout << "Usage:  --gen <Lx> <Ly> <Nx> <Ny> <k3> <k4> | --file \n\t--print | --out (<path>) | --vtk <filename>" << endl;
        cout << "Commands:" << endl;
        cout << "\t--gen                     Generate mesh" << endl; 
        cout << "\t--file                    Read mesh from file(in current directory)" << endl; 
        cout << "\t--print                   Print mesh to stdout" << endl;
        cout << "\t--out                     Print mesh to file" << endl;
        cout << "\t--vtk                     Print mesh in \" .vtk \" format to file" << endl;
        cout << "Options:" << endl; 
        cout << "\t<Lx>                      Mesh X length" << endl;
        cout << "\t<Ly>                      Mesh Y length" << endl;
        cout << "\t<Nx>                      Number of X nodes" << endl;
        cout << "\t<Ny>                      Number of Y nodes" << endl;
        cout << "\t<k3>                      Number of squares in sequence" << endl;
        cout << "\t<k4>                      Number of triangles in sequence" << endl;
        cout << "\t<path>                    Full output path, not necessary. By default, write to current directory." << endl;
        cout << "\t<filename>                Vtk output file name" << endl;
        return 0;
    }
    try {
//creating mesh
        if (!strcmp(argv[1], "--gen")){
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
            topoEN = topos::build_topoEN(Nx, Ny, k3, k4, nE);
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
        else if (!strcmp(argv[1], "--file")){
            if (read_file(C, topoEN,topoSN))
                throw 1;
            type = 1;
            topoNE = topos::build_topoNE(topoEN);
            topoNS = topos::build_topoNS(topoSN);
        }
        else
        {
            cout << "Mesh builder v0.1 beta" << endl;
            cout << "Usage:  --gen <Lx> <Ly> <Nx> <Ny> <k3> <k4> | --file \n\t--print | --out (<path>) | --vtk <filename>" << endl;
            cout << "Commands:" << endl;
            cout << "\t--gen                     Generate mesh" << endl; 
            cout << "\t--file                    Read mesh from file(in current directory)" << endl; 
            cout << "\t--print                   Print mesh to stdout" << endl;
            cout << "\t--out                     Print mesh to file" << endl;
            cout << "\t--vtk                     Print mesh in \" .vtk \" format to file" << endl;
            cout << "Options:" << endl; 
            cout << "\t<Lx>                      Mesh X length" << endl;
            cout << "\t<Ly>                      Mesh Y length" << endl;
            cout << "\t<Nx>                      Number of X nodes" << endl;
            cout << "\t<Ny>                      Number of Y nodes" << endl;
            cout << "\t<k3>                      Number of squares in sequence" << endl;
            cout << "\t<k4>                      Number of triangles in sequence" << endl;
            cout << "\t<path>                    Full output path, not necessary. By default, write to current directory." << endl;
            cout << "\t<filename>                Vtk output file name" << endl;
            return 1;
        }

        cout << "\nMemory:\n" << endl;
        cout << "\tC: " <<  (C).getTotalSize() * sizeof(double) << " bytes" << endl;
        cout << "\ttopoEN: " <<  (topoEN).getTotalSize() * sizeof(int) << " bytes" << endl;
        cout << "\ttopoNE: " <<  (topoNE).getTotalSize() * sizeof(int) << " bytes" << endl;
        cout << "\ttopoSN: " <<  (topoSN).getTotalSize() * sizeof(int) << " bytes" << endl;
        cout << "\ttopoNS: " <<  (topoNS).getTotalSize() * sizeof(int) << " bytes" << endl;


//output mesh
        if (!strcmp(argv[argc-2], "--vtk"))
        {
            char* filename = argv[argc-1];
            strcat(filename,".vtk");
            vtkGenerator<double,int> vtk(filename);
            vtk.printInVTK(nN,C,topoEN);
        }
        else if (!strcmp(argv[argc-2], "--out")){
            write_file(C,topoEN,topoSN,argv[argc-1]);
        }
        else if (!strcmp(argv[argc-1], "--out")){
            write_file(C,topoEN,topoSN);
        }
        else if (!strcmp(argv[argc-1], "--print")){
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
        else
        {
            cout << "Mesh builder v0.1 beta" << endl;
            cout << "Usage:  --gen <Lx> <Ly> <Nx> <Ny> <k3> <k4> | --file \n\t--print | --out (<path>) | --vtk <filename>" << endl;
            cout << "Commands:" << endl;
            cout << "\t--gen                     Generate mesh" << endl; 
            cout << "\t--file                    Read mesh from file(in current directory)" << endl; 
            cout << "\t--print                   Print mesh to stdout" << endl;
            cout << "\t--out                     Print mesh to file" << endl;
            cout << "\t--vtk                     Print mesh in \" .vtk \" format to file" << endl;
            cout << "Options:" << endl; 
            cout << "\t<Lx>                      Mesh X length" << endl;
            cout << "\t<Ly>                      Mesh Y length" << endl;
            cout << "\t<Nx>                      Number of X nodes" << endl;
            cout << "\t<Ny>                      Number of Y nodes" << endl;
            cout << "\t<k3>                      Number of squares in sequence" << endl;
            cout << "\t<k4>                      Number of triangles in sequence" << endl;
            cout << "\t<path>                    Full output path, not necessary. By default, write to current directory." << endl;
            cout << "\t<filename>                Vtk output file name" << endl;
            return 0;
        }
    }
    catch(...){
        cout << "Error occured!" << endl;
        return 1;
    }
    
    return 0;
}