#include <iostream>
#include <fstream>
#include <list>
#include <iterator>
#include <vector>
#include <string.h>
#include <omp.h>
#include <cstring>
#include <cstdlib>

using namespace std;

//Функция количества цифр в числе
int count_digit(int x){
    int i = 0;
    if (x == 0)
        return 1;
    while (x > 0){
        x /= 10;
        i++;
    }
    return i;
}

//Realizing one-dimensional effective-accessible and storing block vector with matrix interface.
template <typename T>
class GridContainer {

    int blockNumber;
    std::vector<int> IA;//offset vector. Size is (N+1)
    std::vector<T> V;//main grid vector
    int countBlockNumber(int vectorSize, const std::vector<int>& blockSizes);
    bool checkSizes(const std::vector<T>& source, const std::vector<int>& blockSizes);

public:

    GridContainer();
    //Creates GridContainer from source vector.
    GridContainer(const std::vector<T>& source, const std::vector<int>& blockSizes);
    GridContainer(GridContainer<T>& source);
    //Can only add elements of vectors which sizes are divided without remainder into {M}
    bool add(const std::vector<T>& extra,const std::vector<int>& blockSizes);
    //Add another GridContainer elements to our container
    bool add(GridContainer<T>& extra);
    GridContainer operator+(GridContainer<T>& extra);
    GridContainer& operator+=(GridContainer<T>& extra);
    T* operator[](int i);
    int getBlocksNumber() const;
    int getBlockSize(int i) const;
    const std::vector<int> &getIa() const;
    void inline printContainer_coord_EN(std::ostream& out = std::cout) {
        for(int i = 0;i<getBlocksNumber();++i) {
            for(int j = 0;j < getBlockSize(i);++j) {
                if (j % 2 == 0) {
                    out << "(";
                    out << operator[](i)[j];
                    out << ", ";
                }
                else {
                    out << operator[](i)[j];
                    out << ")  ";
                }
            }
            out << std::endl;
        }
    }
    void inline printContainer(std::ostream& out = std::cout) {
        for(int i = 0;i<getBlocksNumber();++i) {
            // 15 всего пробелов
            for (int j = 0; j < 7; j++){
                out << " ";
            }
            out << i;
            for (int j = 0; j < 8 - count_digit(i); j++){
                out << " ";
            }
            out << "|    ";
            for(int j = 0;j < getBlockSize(i);++j)
            {
                out << operator[](i)[j];
                if (j != getBlockSize(i) - 1)
                    cout << ",";
            }
            out << std::endl;
        }
    }
};

template <typename T>
int GridContainer<T>::countBlockNumber(int vectorSize, const std::vector<int>& blockSizes) {

    int blockNum = 0;
    int delta = 0;

    for(int i = 0,j = 0;i<vectorSize;++i) {
        delta++;
        if (delta == blockSizes[j]) {
            blockNum++;
            j++;
            delta = 0;
        }
    }

    return blockNum;
}

template <typename T>
bool GridContainer<T>::checkSizes(const std::vector<T>& source, const std::vector<int>& blockSizes) {

    int blSizeElems = 0;
    for(typename std::vector<int>::const_iterator it = blockSizes.begin(); it != blockSizes.end();++it) {
        blSizeElems += *it;
    }

    return (source.size() == blSizeElems);
}

//default constructor
template <typename T>
GridContainer<T>::GridContainer() {
    blockNumber = 0;
}

//Creates GridContainer from source vector.
template <typename T>
GridContainer<T>::GridContainer(const std::vector<T>& source, const std::vector<int>& blockSizes) {
    if (!checkSizes(source,blockSizes)) {
        std::cerr << "Error: Incompatible sizes!" << std::endl;
    }
    else {
        blockNumber = countBlockNumber(source.size(),blockSizes);
        IA.reserve(blockNumber);
        V.reserve(source.size());
        for(typename std::vector<T>::const_iterator it = source.begin(); it != source.end();++it) {
            V.push_back(*it);
        }
        int offset = 0;
        IA.push_back(offset);
        for(typename std::vector<int>::const_iterator it = blockSizes.begin(); it != blockSizes.end();++it) {
            offset += *it;
            IA.push_back(offset);
        }
    }
}

template <typename T>
GridContainer<T>::GridContainer(GridContainer<T>& source){
    blockNumber = source.getBlocksNumber();
    V.resize(source.getIa()[source.getIa().size() - 1]);//allocating
    IA.resize(source.getIa().size());
    for(int i = 0;i<source.getIa().size();++i) {
        IA[i] = source.getIa()[i];
    }
    for(int i = 0;i<source.getBlocksNumber();++i) {

        for(int j = 0;j<source.getBlockSize(i);++j) {
            operator[](i)[j] = source[i][j];
        }
    }
}

//Can only add elements of vectors which sizes are divided without remainder into {M}
template <typename T>
bool GridContainer<T>::add(const std::vector<T>& extra,const std::vector<int>& blockSizes) {

    if (!checkSizes(extra,blockSizes)) {
        std::cerr << "Error: Incompatible sizes!" << std::endl;
        return false;
    }
    else {
        this->blockNumber += countBlockNumber(extra.size(),blockSizes);
        V.reserve(V.size() + extra.size());
        IA.reserve(blockNumber + blockSizes.size());
        for(typename std::vector<T>::const_iterator it = extra.begin(); it != extra.end();++it) {
            V.push_back(*it);
        }
        int offset = IA[IA.size() - 1];
        for(typename std::vector<int>::const_iterator it = blockSizes.begin(); it != blockSizes.end();++it) {
            offset += *it;
            IA.push_back(offset);
        }
        return true;
    }
}

//Add another GridContainer elements to our container
template <typename T>
bool GridContainer<T>::add(GridContainer<T>& extra) {

    int oldBlocksNum = blockNumber;
    int oldOffsetNum = IA.size();
    int offset = IA[IA.size() - 1];
    this->blockNumber += extra.getBlocksNumber();
    IA.resize(IA.size() + extra.getIa().size() - 1);
    V.resize(V.size() + extra.getIa()[extra.getIa().size() - 1]);

    for(int i = oldOffsetNum;i<IA.size();++i) {
        offset += extra.getBlockSize(i - oldOffsetNum );//first if zero
        IA[i] = offset;
    }

    for(int i = oldBlocksNum;i<blockNumber;++i) {
        for(int j = 0;j<getBlockSize(i);++j) {
            operator[](i)[j] = extra[i - oldBlocksNum][j];
        }
    }

    return true;
}

template <typename T>
GridContainer<T> GridContainer<T>::operator+(GridContainer<T>& extra) {
    GridContainer temp_grid(*this);
    temp_grid.add(extra);

    return temp_grid;
}

template <typename T>
GridContainer<T>& GridContainer<T>::operator+=(GridContainer<T>& extra) {
    (*this).add(extra);

    return *this;
}

template <typename T>
T* GridContainer<T>::operator[](int i) {

    if (i >= getBlocksNumber()) {
        std::cerr << "Error: Index out of bounds" << std::endl;
        return nullptr;
    } else {
        return &V[0] + IA[i];
    }
}

template <typename T>
int GridContainer<T>::getBlocksNumber() const {

    return blockNumber;
}

template <typename T>
int GridContainer<T>::getBlockSize(int i) const {
    return IA[i+1] - IA[i];
}

template <typename T>
const std::vector<int>&  GridContainer<T>::getIa() const {
    return IA;
}



template <typename T>
class ConstGridContainer {

    int blocksSize;
    std::vector<T> V;//main grid vector

public:

    ConstGridContainer();
    ConstGridContainer(int num):blocksSize(num){};
    ConstGridContainer(const std::vector<T>& source,int blockSize);
    ConstGridContainer(ConstGridContainer<T>& source);
    void set_blocksize(int num);
    //Can only add elements of vectors which sizes are divided without remainder into {M}
    bool add(const std::vector<T>& extra);
    //Add another GridContainer elements to our container
    bool add(ConstGridContainer<T>& extra);
    ConstGridContainer operator+(ConstGridContainer<T>& extra);
    ConstGridContainer& operator+=(ConstGridContainer<T>& extra);
    T* operator[](int i);
    int getBlocksNumber() const;
    int getBlockSize() const;
    void inline printContainer(std::ostream& out = std::cout) {
        for(int i = 0;i<getBlocksNumber();++i) {
            for(int j = 0;j < blocksSize;++j) {
                if (j % 2 == 0) {
                    out << "(";
                    out << operator[](i)[j];
                    out << ", ";
                }
                else {
                    out << operator[](i)[j];
                    out << ")  ";
                }
            }
            out << std::endl;
        }
    }
};

//default constructor
template <typename T>
ConstGridContainer<T>::ConstGridContainer() {

    blocksSize = 0;
}

template <typename T>
ConstGridContainer<T>::ConstGridContainer(const std::vector<T>& source,int blockSize){

    blocksSize = blockSize;
    V.reserve(source.size());
    for(typename std::vector<T>::const_iterator it = source.begin(); it != source.end();++it) {
        V.push_back(*it);
    }
}

template <typename T>
ConstGridContainer<T>::ConstGridContainer(ConstGridContainer<T>& source){

    blocksSize = source.getBlockSize();
    V.resize(source.getBlockSize()*source.getBlocksNumber());//allocating
    for(int i = 0;i<source.getBlocksNumber();++i) {
        for(int j = 0;j<getBlockSize();++j) {
            operator[](i)[j] = source[i][j];
        }
    }
}

template <typename T>
void ConstGridContainer<T>::set_blocksize(int num){
    if (blocksSize == 0)
        blocksSize = num;
}

//Can only add elements of vectors which sizes are divided without remainder into {M}
template <typename T>
bool ConstGridContainer<T>::add(const std::vector<T>& extra) {

    if (extra.size() % blocksSize != 0) {
        std::cerr << "Error: Vector size must be a multiple of " <<  blocksSize << std::endl;
        return false;
    }
    else {
        V.reserve(extra.size());
        for(typename std::vector<T>::const_iterator it = extra.begin(); it != extra.end();++it) {
            V.push_back(*it);
        }
        return true;
    }
}

//Add another GridContainer elements to our container
template <typename T>
bool ConstGridContainer<T>::add(ConstGridContainer<T>& extra) {

    if (extra.getBlockSize() != getBlockSize()) {
        std::cerr << "Error: Incompatible block sizes!" << std::endl;
        return false;
    }
    else {
        int odlBlocksNum = getBlocksNumber();
        V.resize((getBlocksNumber() + extra.getBlocksNumber()) * getBlockSize());
        for(int i = odlBlocksNum;i<getBlocksNumber();++i) {
            for(int j = 0;j<getBlockSize();++j) {
                operator[](i)[j] = extra[i - odlBlocksNum][j];
            }
        }
        return true;
    }
}

template <typename T>
ConstGridContainer<T> ConstGridContainer<T>::operator+(ConstGridContainer<T>& extra) {

    ConstGridContainer temp_grid(*this);
    temp_grid.add(extra);

    return temp_grid;
}

template <typename T>
ConstGridContainer<T>& ConstGridContainer<T>::operator+=(ConstGridContainer<T>& extra) {

    (*this).add(extra);

    return *this;
}

template <typename T>
T* ConstGridContainer<T>::operator[](int i) {

    if (i >= getBlocksNumber()) {
        std::cerr << "Error: Index out of bounds" << std::endl;
        return nullptr;
    }
    else {
        return &V[0] + getBlockSize()*i;
    }
}

template <typename T>
int ConstGridContainer<T>::getBlocksNumber() const {

    return V.size()/getBlockSize();
}

template <typename T>
int ConstGridContainer<T>::getBlockSize() const {

    return blocksSize;
}


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
GridContainer<int> build_topoNE(GridContainer<int> topoEN){
    vector<int> BlockSize;
    vector<int> temp;
    GridContainer<int> topoNE(temp, BlockSize);
    int nE = topoEN.getBlocksNumber();

    int nN = 0;
    // Находим число точек (max в topoEN + 1)
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

// Построение topoSN
GridContainer<int> build_topoSN(int Nx, int Ny){
    vector<int> BlockSize;
    vector<int> temp;
    GridContainer<int> topoSN(temp, BlockSize);

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

// Построение topoNS
GridContainer<int> build_topoNS(int Nx, int Ny){
    vector<int> BlockSize;
    vector<int> temp;
    GridContainer<int> topoNS(temp, BlockSize);

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

// Запись в файл сеточных данных
void write_file(ConstGridContainer<double> C, GridContainer<int> topoEN, GridContainer<int> topoSN){
    ofstream fout;

    fout.open("mesh.txt");
    fout << "Nn " << C.getBlocksNumber() << '\n';
    fout << "Nt " << topoEN.getBlocksNumber() << '\n';
    fout << "NFaceBC " << topoSN.getBlocksNumber() << '\n';
    fout << "NumCoords " << C.getBlockSize() << '\n';
    fout.close();

    fout.open("coordinate.msh");
    for (int i = 0; i < C.getBlocksNumber(); i++)
        fout << C[i][0] << " " << C[i][1] << '\n';
    fout.close();

    fout.open("topo.msh");
    for (int i = 0; i < topoEN.getBlocksNumber(); i++) {
        fout << topoEN.getBlockSize(i);
        for (int j = 0; j < topoEN.getBlockSize(i); j++)
            fout << " " << topoEN[i][j];
        fout << '\n';
    }
    fout.close();

    fout.open("bctopo.msh");
    for (int i = 0; i < topoSN.getBlocksNumber(); i++){
        fout << topoSN.getBlockSize(i) - 1;
        for (int j = 1; j < topoSN.getBlockSize(i); j++)
            fout << " " << topoSN[i][j];
        fout << " " << topoSN[i][0];
        fout << '\n';
    }
    fout.close();
}

// Чтение из файла
int read_file(ConstGridContainer<double> &C, GridContainer<int> &topoEN, GridContainer<int> &topoSN){
    ifstream fin;
    int nN, nE, NFaceBC, NumCoords, count_node;
    char str;
    vector<int> BlockSize;


    fin.open("mesh.txt");
    if (fin.is_open()){
        for (int i = 0; i < 2; i++) fin >> str;
        fin >> nN;
        for (int i = 0; i < 2; i++) fin >> str;
        fin >> nE;
        for (int i = 0; i < 7; i++) fin >> str;
        fin >> NFaceBC;
        for (int i = 0; i < 9; i++) fin >> str;
        fin >> NumCoords;
    }
    else {
        cout << '\n' << "File \"mesh.txt\" not found" << '\n';
        return 1;
    }
    fin.close();


    fin.open("coordinate.msh");
    C.set_blocksize(NumCoords);
    if (fin.is_open()){
        vector<double> temp;
        double coord;

        for (int i = 0; i < nN; i++){
            for (int j = 0; j < NumCoords; j++){
                fin >> coord;
                temp.push_back(coord);
            }
            C.add(temp);
            temp.clear();
        }
    }
    else {
        cout << '\n' << "File \"coordinate.msh\" not found" << '\n';
        return 1;
    }
    fin.close();

    fin.open("topo.msh");
    if (fin.is_open()){
        vector<int> temp;
        BlockSize.push_back(0);
        topoEN.add(temp, BlockSize);
        BlockSize.clear();
        int coord;

        for (int i = 0; i < nE; i++){
            fin >> count_node;
            BlockSize.push_back(count_node);

            for (int j = 0; j < count_node; j++){
                fin >> coord;
                temp.push_back(coord);
            };

            topoEN.add(temp, BlockSize);

            temp.clear();
            BlockSize.clear();
        }
    }
    else {
        cout << '\n' << "File \"topo.msh\" not found" << '\n';
        return 1;
    }
    fin.close();

    fin.open("bctopo.msh");
    if (fin.is_open()){
        vector<int> temp;
        BlockSize.push_back(0);
        topoSN.add(temp, BlockSize);
        BlockSize.clear();
        int number_edge;
        int coord;

        for (int i = 0; i < NFaceBC; i++){
            fin >> count_node;
            BlockSize.push_back(count_node+1);

            temp.push_back(0);
            for (int j = 0; j < count_node; j++){
                fin >> coord;
                temp.push_back(coord);
            };
            fin >> number_edge;
            temp[0] = number_edge;

            topoSN.add(temp, BlockSize);

            temp.clear();
            BlockSize.clear();
        }

    }
    else {
        cout << '\n' << "File \"bctopo.msh\" not found" << '\n';
        return 1;
    }
    fin.close();


    return 0;
}

int main(int argc, char **argv) {
    int Nx, Ny, k3, k4, nN, nE;
    double Lx, Ly;
    ConstGridContainer<double> C;
    GridContainer<int> topoEN;
    GridContainer<int> topoNE;
    GridContainer<int> topoSN;
    GridContainer<int> topoNS;
    int type = 0;

    double start, end;

    if ((argc == 2) && (!strcmp(argv[1], "-help"))) {
        cout << "Usage:\n-gen Lx Ly Nx Ny k3 k4\n\t or\n-file (to read data from files)" << endl;
        return 1;
    }
    else if ((argc == 2) && (!strcmp(argv[1], "-file"))){
        if (read_file(C, topoEN,topoSN))
            return 1;
        type = 1;
        topoNE = build_topoNE(topoEN);
    }
    else if ((argc == 8) && (!strcmp(argv[1], "-gen"))){
        Lx = atoi(argv[2]);
        Ly = atoi(argv[3]);
        Nx = atoi(argv[4]);
        Ny = atoi(argv[5]);
        k3 = atoi(argv[6]);
        k4 = atoi(argv[7]);
        if (k3 * k3 + k4 * k4 == 0){
            cout << "Incorrect data ((k3 or k4) != 0)\n";
            return 1;
        }

        C.set_blocksize(2);
        nN = Nx *  Ny;
        nE = num_elem(Nx, Ny, k3, k4);

        start = omp_get_wtime();
        build_coord(C, Lx, Ly, Nx, Ny);
        end = omp_get_wtime();
        cout << "build C: " << end - start << "sec" << endl;

        start = omp_get_wtime();
        topoEN = build_topoEN(C, Nx, Ny, k3, k4, nE);
        end = omp_get_wtime();
        cout << "build topoEN: " << end - start << "sec" << endl;

        start = omp_get_wtime();
        topoNE = build_topoNE(topoEN);
        end = omp_get_wtime();
        cout << "build topoNE: " << end - start << "sec" << endl;

        start = omp_get_wtime();
        topoSN = build_topoSN(Nx,Ny);
        end = omp_get_wtime();
        cout << "build topoSN: " << end - start << "sec" << endl;

        start = omp_get_wtime();
        topoNS = build_topoNS(Nx,Ny);
        end = omp_get_wtime();
        cout << "build topoSN: " << end - start << "sec" << endl;
    }
    else {
        cout << "Use -help" << endl;
        return 1;
    }

    string req, toponame;
    cout << "Аbilities:\n\tprint <toponame>\n\toutput_mesh\n\toutput (to write data in files)\n\tq - quit\n";

    while(true){
        cin >> req;

        if (req == "output_mesh"){
            if (!type) {
                cout << "\nGetting mesh: ";
                draw_grid(Nx - 1, Ny - 1, k3, k4);
            }
            else
                cout << "\nUnfortunately, mesh doesn't print in this mode\n\n";
        }
        else if (req == "output"){
            write_file(C,topoEN,topoSN);
            cout << "\nDone\n\n";
        }
        else if (req == "print"){
            cin >> toponame;
            if (toponame == "topoEN") {
                topoEN.printContainer();
                cout << '\n';
            }
            if (toponame == "topoNE") {
                topoNE.printContainer();
                cout << '\n';
            }
            if (toponame == "C") {
                C.printContainer();
                cout << '\n';
            }
            if (toponame == "topoSN") {
                topoSN.printContainer();
                cout << '\n';
            }
            //if (toponame == "topoNS")
        }
        else if (req == "q") return 0;
    }

    /*
    C.set_blocksize(2);

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
        if (!strcmp(argv[1],"file")){

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
    start = omp_get_wtime();
    build_coord(C, Lx, Ly, Nx, Ny);
    end = omp_get_wtime();
    C.printContainer();
    cout << "build C: " << end - start << "sec" << endl;
    cout << "\nПолученная сетка: ";
    draw_grid(Nx - 1, Ny - 1, k3, k4);

    cout << "topoEN: \n";
    cout << "Номер элемента | Номера точек \n";
    start = omp_get_wtime();
    topoEN = build_topoEN(C, Nx, Ny, k3, k4, nE);
    end = omp_get_wtime();
    topoEN.printContainer();
    cout << "build topoEN: " << end - start << "sec" << endl;

    cout << "\ntopoNE: \n";
    cout << "Номер точки    | Номера элементов \n";
    start = omp_get_wtime();
    topoNE = build_topoNE(topoEN);
    end = omp_get_wtime();
    topoNE.printContainer();
    cout << "build topoNE: " << end - start << "sec" << endl;

    cout << "\ntopoSN: \n";
    cout << "Номер сегмента | Номер поверхности, номера точек \n";
    start = omp_get_wtime();
    topoSN = build_topoSN(Nx,Ny);
    end = omp_get_wtime();
    topoSN.printContainer();
    cout << "build topoSN: " << end - start << "sec" << endl;

    cout << "\ntopoNS: \n";
    cout << "Номер точки    | Номера сегментов \n";
    start = omp_get_wtime();
    topoNS = build_topoNS(Nx,Ny);
    end = omp_get_wtime();
    topoNS.printContainer();
    cout << "build topoSN: " << end - start << "sec" << endl;

    cout << "\nКоличество точек: " << nN << '\n';
    cout << "Количество элементов: " << nE;
    cout << "\n";

    write_file(C,topoEN, topoSN);


    return 0;
    */
}
