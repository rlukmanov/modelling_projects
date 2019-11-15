#include "FixedSizeMeshContainer.cpp"
#include "VariableSizeMeshContainer.cpp"
#include <omp.h>

#pragma once

namespace topos
{

    //coords
    void build_coord(FixedSizeMeshContainer<double>& C, double Lx, double Ly, int Nx, int Ny){
        vector<double> temp;

        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                temp.push_back((Lx / (Nx - 1)) * j);
                temp.push_back((Ly / (Ny - 1)) * i);
                C.add(temp);
                temp.clear();
            }
        }
    }
    
    // topoEN
    VariableSizeMeshContainer<int> build_topoEN(int Nx, int Ny, int k3, int k4, int nE){
        vector<int> BlockSize;
        BlockSize.push_back(0);
        vector<int> temp;
        VariableSizeMeshContainer<int> topoEN;
        bool figure = k3 > 0 ? false : true; // 0 - triangle, 1 - square
        int count_figure = figure ? k4 : k3;
        int EN_i = 0;
        int temp_nE;

        temp_nE = nE;
        topoEN.add(temp, BlockSize);
        BlockSize.clear();

        while (temp_nE>0) {
            if (!figure) {
                temp.push_back(EN_i);
                temp.push_back(EN_i+1);
                temp.push_back(EN_i+Nx);

                BlockSize.push_back(3);


                temp.push_back(EN_i+1);
                temp.push_back(EN_i+1+Nx);
                temp.push_back(EN_i+Nx);

                BlockSize.push_back(3);


                temp_nE-=2;
            }
            else 
            {
                temp.push_back(EN_i);
                temp.push_back(EN_i+1);
                temp.push_back(EN_i+1+Nx);
                temp.push_back(EN_i+Nx);
                
                BlockSize.push_back(4);

                temp_nE--;
            }
            count_figure--;
            if (count_figure == 0) {
                figure = !figure;
                count_figure = figure ? k4 : k3;
                if (count_figure == 0) {
                    figure = !figure;
                    count_figure = figure ? k4 : k3;
                }
            }
            EN_i++;
            if (EN_i % Nx == Nx - 1)
                EN_i++;
        }

        topoEN.add(temp, BlockSize);
        temp.clear();
        BlockSize.clear();

        return topoEN;
    }

    // topoNE
    VariableSizeMeshContainer<int> build_topoNE(const VariableSizeMeshContainer<int>& topoEN){
        vector<int> BlockSize;
        vector<int> temp;
        VariableSizeMeshContainer<int> topoNE(temp, BlockSize);
        int nE = topoEN.getBlockNumber();

        int nN = 0;

        for (int i = 0; i < nE; i++){
            for (int j = 0; j < topoEN.getBlockSize(i); j++)
                if (topoEN[i][j]>nN)
                    nN = topoEN[i][j];
        }
        nN++;

        int count_mas[nN];
        int count_i_mas[nN];

        for (int i = 0; i < nN; i++){
            count_mas[i] = 0;
            count_i_mas[i] = 0;
        }


        for (int i = 0; i < nE; i++){
            for (int j = 0; j < topoEN.getBlockSize(i); j++) {
                count_mas[topoEN[i][j]]++;
                count_i_mas[topoEN[i][j]]++;
            }
        }

        for (int i = 0; i < nN; i++){
            BlockSize.push_back(count_mas[i]);
            for (int j = 0; j < count_mas[i]; j++)
                temp.push_back(0);
        }

        topoNE.add(temp,BlockSize);
        BlockSize.clear();
        temp.clear();

        for (int i = 0; i < nE; i++){
            for (int j = 0; j < topoEN.getBlockSize(i); j++){
                topoNE[topoEN[i][j]][count_mas[topoEN[i][j]] - count_i_mas[topoEN[i][j]]] = i;
                count_i_mas[topoEN[i][j]]--;
            }
        }

        return topoNE;
    }

    // topoSN
    VariableSizeMeshContainer<int> build_topoSN(int Nx, int Ny){
        vector<int> BlockSize;
        vector<int> temp;
        VariableSizeMeshContainer<int> topoSN(temp, BlockSize);

        for(int i = 0; i < Nx - 1; ++i) {
            temp.push_back(0);
            temp.push_back(i);
            temp.push_back(i + 1);
            BlockSize.push_back(3);

            topoSN.add(temp, BlockSize);

            temp.clear();
            BlockSize.clear();
        }

        for(int i = 0; i < Ny - 1; ++i) {
            temp.push_back(1);
            temp.push_back(Nx - 1 + i);
            temp.push_back(Nx - 1 + i + 1);
            BlockSize.push_back(3);

            topoSN.add(temp, BlockSize);

            temp.clear();
            BlockSize.clear();
        }

        for(int i = 0; i < Nx - 1; ++i) {
            temp.push_back(2);
            temp.push_back(Nx + Ny - 2 + i);
            temp.push_back(Nx + Ny - 2 + i + 1);
            BlockSize.push_back(3);

            topoSN.add(temp, BlockSize);

            temp.clear();
            BlockSize.clear();
        }

        for(int i = 0; i < Ny - 1; ++i) {
            temp.push_back(3);
            temp.push_back(Nx + Nx + Ny - 3 + i);

            if (i + 1 == Ny - 1)
                temp.push_back(0);
            else
                temp.push_back(Nx + Nx + Ny - 3 + i + 1);
            BlockSize.push_back(3);

            topoSN.add(temp, BlockSize);

            temp.clear();
            BlockSize.clear();
        }

        return topoSN;
    }

    // topoNS
    VariableSizeMeshContainer<int> build_topoNS(const VariableSizeMeshContainer<int>& topo){
        vector<int> BlockSize;
        vector<int> temp;
        int Nx, Ny;
        VariableSizeMeshContainer<int> topoNS(temp, BlockSize);

        for(Nx = 1; topo[Nx - 1][0] == 0; ++Nx) {}

        for(Ny = 1; topo[Nx + Ny - 2][0] == 1; ++Ny) {}


        for(int i = 0; i < (Nx + Ny - 2) * 2; ++i) {
            if (i == 0) 
                temp.push_back(2*(Nx + Ny - 2) - 1);
            else 
                temp.push_back(i - 1);

            temp.push_back(i);
            BlockSize.push_back(2);

            topoNS.add(temp, BlockSize);

            temp.clear();        
            BlockSize.clear();
        }

        return topoNS;
    }
}