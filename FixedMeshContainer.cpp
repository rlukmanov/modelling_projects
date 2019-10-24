#include <vector>
#include <iostream>
#ifndef __CONST_GRID
#define __CONST_GRID

//Realizing one-dimensional effective-accessible and storing block vector with matrix interface.
template <typename T>
class FixedSizeMeshContainer
{

    int blocksSize;
    std::vector<T> V;//main grid vector

public:

    FixedSizeMeshContainer();

    FixedSizeMeshContainer(int num);

    FixedSizeMeshContainer(const std::vector<T>& source,int blockSize);

    FixedSizeMeshContainer(FixedSizeMeshContainer<T>& source);

    //Can only add elements of vectors which sizes are divided without remainder into {M}
    bool add(const std::vector<T>& extra);

    //Add another GridContainer elements to our container
    bool add(FixedSizeMeshContainer<T>& extra);


    FixedSizeMeshContainer operator+(FixedSizeMeshContainer<T>& extra);

    FixedSizeMeshContainer& operator+=(FixedSizeMeshContainer<T>& extra);

    T* operator[](int i);

    int getBlocksNumber() const;

    int getBlockSize() const;

    void inline printContainer(std::ostream& out = std::cout)
    {
        for(int i = 0;i<getBlocksNumber();++i) 
        {
            for(int j = 0;j < blocksSize;++j)
            {
                if (j % 2 == 0)
                {
                    out << "(";
                    out << operator[](i)[j];
                    out << ", ";
                }
                else
                {
                    out << operator[](i)[j];
                    out << ")  ";
                }
            }
            out << std::endl;
        }
    }

    // void inline printContainer(std::ostream& out = std::cout)
    // {
    //     for(int i = 0;i<getBlocksNumber();++i)
    //     {
    //         for(int j = 0;j < blocksSize;++j)
    //         {
    //             out << operator[](i)[j];
    //             out << " ";
    //         }
    //         out << std::endl;
    //     }
    // }
};

//default constructor
template <typename T>
FixedSizeMeshContainer<T>::FixedSizeMeshContainer()
{
    blocksSize = 0;
}

template<typename T>
FixedSizeMeshContainer<T>::FixedSizeMeshContainer(int num):blocksSize(num){};

template <typename T>
FixedSizeMeshContainer<T>::FixedSizeMeshContainer(const std::vector<T>& source,int blockSize)
{
    blocksSize = blockSize;
    V.reserve(source.size());
    for(typename std::vector<T>::const_iterator it = source.begin(); it != source.end();++it)
    {
        V.push_back(*it);
    }
}

template <typename T>
FixedSizeMeshContainer<T>::FixedSizeMeshContainer(FixedSizeMeshContainer<T>& source)
{
    blocksSize = source.getBlockSize();
    V.resize(source.getBlockSize()*source.getBlocksNumber());//allocating

    V = source.V;
}


//Can only add elements of vectors which sizes are divided without remainder into {M}
template <typename T>
bool FixedSizeMeshContainer<T>::add(const std::vector<T>& extra)
{
    if (extra.size() % blocksSize != 0)
    {
        std::cerr << "Error: Vector size must be a multiple of " <<  blocksSize << std::endl;
        return false;
    } 
    else
    {
        V.reserve(extra.size());
        for(typename std::vector<T>::const_iterator it = extra.begin(); it != extra.end();++it)  
        {
            V.push_back(*it);
        }
        return true;
    }
}

//Add another GridContainer elements to our container
template <typename T>
bool FixedSizeMeshContainer<T>::add(FixedSizeMeshContainer<T>& extra)
{
    if (extra.getBlockSize() != getBlockSize())
    {
        std::cerr << "Error: Incompatible block sizes!" << std::endl;
        return false;
    }
    else
    {
        int odlBlocksNum = getBlocksNumber();
        V.resize((getBlocksNumber() + extra.getBlocksNumber()) * getBlockSize());
        for(int i = odlBlocksNum;i<getBlocksNumber();++i)
        {
            for(int j = 0;j<getBlockSize();++j)
            {
                operator[](i)[j] = extra[i - odlBlocksNum][j];
            }
        }
        return true;
    }
}

template <typename T>
FixedSizeMeshContainer<T> FixedSizeMeshContainer<T>::operator+(FixedSizeMeshContainer<T>& extra)
{
    FixedSizeMeshContainer temp_grid(*this);
    temp_grid.add(extra);
    return temp_grid;
}

template <typename T>
FixedSizeMeshContainer<T>& FixedSizeMeshContainer<T>::operator+=(FixedSizeMeshContainer<T>& extra)
{
    (*this).add(extra);
    return *this;
}

template <typename T>
T* FixedSizeMeshContainer<T>::operator[](int i)
{
    if (i >= getBlocksNumber())
    {
        std::cerr << "Error: Index out of bounds" << std::endl;
        return nullptr;
    } 
    else
    {
        return &V[0] + getBlockSize()*i;
    }
}

template <typename T>
int FixedSizeMeshContainer<T>::getBlocksNumber() const
{
    return V.size()/getBlockSize();
}

template <typename T>
int FixedSizeMeshContainer<T>::getBlockSize() const
{
    return blocksSize;
}


#endif