#include <iostream>
#include <list>
#include <iterator>
#include <vector>
#include <cstdlib>
#include <cstring>
#include "GridContainer.cpp"
#include "ConstGridContainer.cpp"

using namespace std;

// Отрисовка сетки
void draw_grid(int Nx, int Ny, int k3, int k4){
    int cur_number_type = k3 > 0 ? k3:k4; // Текущее количество элементов
    int type_elem = k3 > 0 ? 0:1; // Тип элемента. 0 - треугольник, 1 квадрат
    int prev_cur_number = cur_number_type;
    int prev_type = type_elem;
    char symb = ' ';

    cout << "\n";

    for (int i = 0; i < Nx; i++) {
        cout << " ";
        for (int j = 0; j < 8; j++)
            cout << "_";
    }
    cout << "\n";

    for (int cur_line = 0; cur_line < Ny; cur_line++) {
        symb = ' ';

        for (int cur_str = 0; cur_str < 4; cur_str++) {
            cout << "|";
            if (cur_str == 3)
                symb = '_';

            type_elem = prev_type;
            cur_number_type = type_elem == 0 ? k3:k4;
            if ((cur_line != 0) && (prev_cur_number != 0 )){
                cur_number_type = prev_cur_number;
            }

            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < 7 - 2 * cur_str; j++)
                    cout << symb;
                if (type_elem == 0) {
                    cout << "̸"; //U+0338
                }
                else {
                    cout << symb;
                }
                for (int j = 0; j < 2 * cur_str; j++)
                    cout << symb;
                cout << "|";

                cur_number_type--;
                if (cur_number_type == 0){
                    type_elem = (type_elem + 1) % 2;
                    cur_number_type = type_elem == 0 ? k3:k4;
                    if (cur_number_type == 0) {
                        type_elem = (type_elem + 1) % 2;
                        cur_number_type = type_elem == 0 ? k3:k4;
                    }
                }
            }
            cout << "\n";
        }
        prev_cur_number = cur_number_type;
        prev_type = type_elem;
    }

    cout << "\n";
}

// Подсчет количества элементов
int num_elem(int Nx, int Ny, int k3, int k4){
    int nE;

    nE = (2 * k3 + k4) * int(((Nx - 1) * (Ny - 1)) / (k3 + k4));
    if ((((Nx - 1) * (Ny - 1)) % (k3 + k4)) >= k3)
        nE += 2 * k3 + (((Nx - 1) * (Ny - 1)) % (k3 + k4)) - k3;
    else
        nE += 2 * (((Nx - 1) * (Ny - 1)) % (k3 + k4));

    return nE;
}

// Построение массива координат
void build_coord(ConstGridContainer<double>& C, double Lx, double Ly, int Nx, int Ny){
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

// Построение topoEN
GridContainer<int> build_topoEN(ConstGridContainer<double> C, int Nx, int Ny, int k3, int k4, int nE){
    vector<int> BlockSize;
    BlockSize.push_back(0);
    vector<int> temp;
    GridContainer<int> topoEN;
    bool figure = k3 > 0 ? false : true; // 0 - тругольник, 1 - квадрат
    int count_figure = figure ? k4 : k3;
    int EN_i = 0;
    int temp_nE;

    temp_nE = nE;
    topoEN.add(temp, BlockSize);
    BlockSize.clear();

    while (temp_nE>0) {
        if (!figure) {
            //temp.push_back(C[EN_i][0]);
            //temp.push_back(C[EN_i][1]);

            //temp.push_back(C[EN_i + 1][0]);
            //temp.push_back(C[EN_i + 1][1]);

            //temp.push_back(C[EN_i + Nx][0]);
            //temp.push_back(C[EN_i + Nx][1]);

            temp.push_back(EN_i);
            temp.push_back(EN_i+1);
            temp.push_back(EN_i+Nx);
            BlockSize.push_back(3);


            //BlockSize.push_back(6);

            topoEN.add(temp, BlockSize);
            temp.clear();
            BlockSize.clear();


            //temp.push_back(C[EN_i + 1][0]);
            //temp.push_back(C[EN_i + 1][1]);

            //temp.push_back(C[EN_i + 1 + Nx][0]);
            //temp.push_back(C[EN_i + 1 + Nx][1]);

            //temp.push_back(C[EN_i + Nx][0]);
            //temp.push_back(C[EN_i + Nx][1]);

            //BlockSize.push_back(6);

            temp.push_back(EN_i+1);
            temp.push_back(EN_i+1+Nx);
            temp.push_back(EN_i+Nx);
            BlockSize.push_back(3);

            topoEN.add(temp, BlockSize);
            temp.clear();
            BlockSize.clear();
            temp_nE-=2;
        } else {
            //temp.push_back(C[EN_i][0]);
            //temp.push_back(C[EN_i][1]);

            //temp.push_back(C[EN_i + 1][0]);
            //temp.push_back(C[EN_i + 1][1]);

            //temp.push_back(C[EN_i + 1 + Nx][0]);
            //temp.push_back(C[EN_i + 1 + Nx][1]);

            //temp.push_back(C[EN_i + Nx][0]);
            //temp.push_back(C[EN_i + Nx][1]);

            //BlockSize.push_back(8);

            temp.push_back(EN_i);
            temp.push_back(EN_i+1);
            temp.push_back(EN_i+1+Nx);
            temp.push_back(EN_i+Nx);
            BlockSize.push_back(4);

            topoEN.add(temp, BlockSize);
            temp.clear();
            BlockSize.clear();
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

    return topoEN;
}

// Построение topoNE
GridContainer<int> build_topoNE(ConstGridContainer<double> C, GridContainer<int> topoEN, int Nx, int Ny, int k3, int k4, int nE){
    vector<int> BlockSize;
    vector<int> temp;
    GridContainer<int> topoNE(temp, BlockSize);
    int nN = Nx * Ny;
    int count_mas[nN];
    int count_i_mas[nN];

    for (int i = 0; i < nN; i++){
        count_mas[i] = 0;
        count_i_mas[i] = 0;
    }

    //Проход для подсчета количества встреч каждой точки
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
        topoNE.add(temp,BlockSize);
        BlockSize.clear();
        temp.clear();
    }

    for (int i = 0; i < nE; i++){
        for (int j = 0; j < topoEN.getBlockSize(i); j++){
            topoNE[topoEN[i][j]][count_mas[topoEN[i][j]] - count_i_mas[topoEN[i][j]]] = i;
            count_i_mas[topoEN[i][j]]--;
        }
    }

    return topoNE;
}

int main(int argc, char **argv) {
    int Nx, Ny, k3, k4, nN, nE;
    double Lx, Ly;
    ConstGridContainer<double> C(2);
    GridContainer<int> topoEN;
    GridContainer<int> topoNE;

    if (argc == 2){
        if (!strcmp(argv[1],"-help")) {
            cout << "Введите длину сетки (Lx),\n";
            cout << "Введите ширину сетки (Ly),\n";
            cout << "Введите количество точек (Nx),\n";
            cout << "Введите количество точек (Ny),\n";
            cout << "Введите количество элементов с треугольным элементом (идут подряд) (k3),\n";
            cout << "Введите количество элементов с квадратным элементом (идут подряд) (k4).\n";
            return 0;
        }
    }
    else if (argc == 7) {
        Lx = atoi(argv[1]);
        Ly = atoi(argv[2]);
        Nx = atoi(argv[3]);
        Ny = atoi(argv[4]);
        k3 = atoi(argv[5]);
        k4 = atoi(argv[6]);
        if (k3 * k3 + k4 * k4 == 0){
            cout << "Введите ненулевое количество треугольников и квадратов\n";
            return 1;
        }
    }
    else {
        cout << "Недостаточное число аргументов. Попробуйте -help\n";
        return 0;
    }

    nN = Nx * Ny;
    nE = num_elem(Nx, Ny, k3, k4);

    cout << "\nМассив координат: \n";
    cout << "x | y\n";
    build_coord(C, Lx, Ly, Nx, Ny);
    C.printContainer();

    cout << "\nПолученная сетка: ";
    draw_grid(Nx - 1, Ny - 1, k3, k4);

    cout << "topoEN: \n";
    cout << "Номер элемента | Номера точек \n";
    topoEN = build_topoEN(C, Nx, Ny, k3, k4, nE);
    topoEN.printContainer();

    cout << "\ntopoNE: \n";
    cout << "Номер точки    | Номера элементов \n";
    topoNE = build_topoNE(C,topoEN,Nx,Ny,k3,k4,nE);
    topoNE.printContainer();

    cout << "\nКоличество точек: " << nN << '\n';
    cout << "Количество элементов: " << nE;
    cout << "\n";

    return 0;
}
