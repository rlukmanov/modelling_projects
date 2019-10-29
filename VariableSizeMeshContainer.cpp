#include <vector>
#include <iostream>
#pragma once

//Realizing one-dimensional effective-accessible and storing block vector with matrix interface.
template <typename T>
class VariableSizeMeshContainer
{
    int blockNumber;
    std::vector<int> IA;//offset vector. Size is (N+1)
    std::vector<T> V;//main grid vector

    int countBlockNumber(int vectorSize, const std::vector<int>& blockSizes);

    bool checkSizes(const std::vector<T>& source, const std::vector<int>& blockSizes);

    int count_digit(int x);

public:

    VariableSizeMeshContainer();
    
    //Creates VariableSizeMeshContainer from source vector.
    VariableSizeMeshContainer(const std::vector<T>& source, const std::vector<int>& blockSizes);

    VariableSizeMeshContainer(const VariableSizeMeshContainer<T>& source);

    //Can only add elements of vectors which sizes are divided without remainder into {M}
    bool add(const std::vector<T>& extra,const std::vector<int>& blockSizes);


    //Add another VariableSizeMeshContainer elements to our container
    bool add(const VariableSizeMeshContainer<T>& extra);

    VariableSizeMeshContainer& operator=(const VariableSizeMeshContainer<T>& target);

    VariableSizeMeshContainer operator+(const VariableSizeMeshContainer<T>& extra);

    VariableSizeMeshContainer& operator+=(const VariableSizeMeshContainer<T>& extra);

    T* operator[](int i);

    const T* operator[](int i) const;

    int getBlockNumber() const;

    int getBlockSize(int i) const;

    // void inline printContainer(std::ostream& out = std::cout)
    // {
    //     for(int i = 0;i<BlocksNumber;++i)
    //     {
    //         for(int j = 0;j < getBlockSize(i);++j)
    //         {
    //             out << operator[](i)[j];
    //             out << " ";
    //         }
    //         out << std::endl;
    //     }
    // }

    void inline printContainer_coord_EN(std::ostream& out = std::cout)
    {
        for(int i = 0;i<blockNumber;++i)
        {
            for(int j = 0;j < getBlockSize(i);++j)
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
    void inline printContainer(std::ostream& out = std::cout)
    {
        for(int i = 0;i<blockNumber;++i)
        {
            // 15 всего пробелов
            for (int j = 0; j < 7; j++)
            {
                out << " ";
            }
            out << i;
            for (int j = 0; j < 8 - count_digit(i); j++)
            {
                out << " ";
            }
            out << "|    ";
            for(int j = 0;j < getBlockSize(i);++j)
            {
                out << operator[](i)[j];
                if (j != getBlockSize(i) - 1)   out << ",";
            }
            out << std::endl;
        }
    }
};

template<typename T>
int VariableSizeMeshContainer<T>::count_digit(int x)
{
    int i = 0;
    if (x == 0) return 1;
    while (x > 0)
    {
        x /= 10;
        i++;
    }
    return i;
}

template <typename T>
int VariableSizeMeshContainer<T>::countBlockNumber(int vectorSize, const std::vector<int>& blockSizes)
{
    int blockNum = 0;
    int delta = 0;
    for(int i = 0,j = 0;i<vectorSize;++i)
    {
        delta++;
        if (delta == blockSizes[j])
        {
            blockNum++;
            j++;
            delta = 0;
        }
    }
    return blockNum;
}

template <typename T>
bool VariableSizeMeshContainer<T>::checkSizes(const std::vector<T>& source, const std::vector<int>& blockSizes)
{
    unsigned int blSizeElems = 0;
    for(typename std::vector<int>::const_iterator it = blockSizes.begin(); it != blockSizes.end();++it)
    {
        blSizeElems += *it;
    }
    return (source.size() == blSizeElems);
}

//default constructor
template <typename T>
VariableSizeMeshContainer<T>::VariableSizeMeshContainer()
{
    blockNumber = 0;
}


//Creates VariableSizeMeshContainer from source vector.
template <typename T>
VariableSizeMeshContainer<T>::VariableSizeMeshContainer(const std::vector<T>& source, const std::vector<int>& blockSizes)//number is size of each vector's block
{
    if (!checkSizes(source,blockSizes))
    {
        std::cerr << "Error: Incompatible sizes!" << std::endl;
    }   
    else
    {
        blockNumber = countBlockNumber(source.size(),blockSizes);
        IA.reserve(blockNumber);
        V.reserve(source.size());
        for(typename std::vector<T>::const_iterator it = source.begin(); it != source.end();++it)
        {
            V.push_back(*it);
        }
        int offset = 0;
        IA.push_back(offset);
        for(typename std::vector<T>::const_iterator it = blockSizes.begin(); it != blockSizes.end();++it)
        {
            offset += *it;
            IA.push_back(offset);
        }
    }
}

template <typename T>
VariableSizeMeshContainer<T>::VariableSizeMeshContainer(const VariableSizeMeshContainer<T>& source)//number is size of each vector's block
{
    blockNumber = source.blockNumber;

    V.resize(source.IA[source.IA.size() - 1]);//allocating
    IA.resize(source.IA.size());

    IA = source.IA;
    V = source.V;
}

//Can only add elements of vectors which sizes are divided without remainder into {M}
template <typename T>
bool VariableSizeMeshContainer<T>::add(const std::vector<T>& extra,const std::vector<int>& blockSizes)
{
    if (!checkSizes(extra,blockSizes))
    {
        std::cerr << "Error: Incompatible sizes!" << std::endl;
        return false;
    } 
    else
    {
        this->blockNumber += countBlockNumber(extra.size(),blockSizes);
        V.reserve(V.size() + extra.size());
        IA.reserve(blockNumber + blockSizes.size());
        for(typename std::vector<T>::const_iterator it = extra.begin(); it != extra.end();++it)
        {
            V.push_back(*it);
        }
        int offset = IA[IA.size() - 1];
        for(typename std::vector<T>::const_iterator it = blockSizes.begin(); it != blockSizes.end();++it)
        {
            offset += *it;
            IA.push_back(offset);
        }
        return true;
    }
}


//Add another VariableSizeMeshContainer elements to our container
template <typename T>
bool VariableSizeMeshContainer<T>::add(const VariableSizeMeshContainer<T>& extra)
{
    int oldBlocksNum = blockNumber;
    int oldOffsetNum = IA.size();
    int offset = IA[IA.size() - 1];

    this->blockNumber += extra.BlocksNumber;
    IA.resize(IA.size() + extra.IA.size() - 1);
    V.resize(V.size() + extra.IA[extra.IA.size() - 1]);


    for(int i = oldOffsetNum;i<IA.size();++i)
    {
        offset += extra.getBlockSize(i - oldOffsetNum );//first if zero
        IA[i] = offset;
    }

    for(int i = oldBlocksNum;i<blockNumber;++i)
    {
        for(int j = 0;j<getBlockSize(i);++j)
        {
            operator[](i)[j] = extra[i - oldBlocksNum][j];
        }
    }
    return true;
}


template <typename T>
VariableSizeMeshContainer<T>& VariableSizeMeshContainer<T>::operator=(const VariableSizeMeshContainer<T>& target)
{
    blockNumber = target.blockNumber;
    
    V.clear();
    IA.clear();

    V.resize(target.IA[target.IA.size() - 1]);//allocating
    IA.resize(target.IA.size());

    IA = target.IA;
    V = target.V;
    
    return (*this);
}

template <typename T>
VariableSizeMeshContainer<T> VariableSizeMeshContainer<T>::operator+(const VariableSizeMeshContainer<T>& extra)
{
    VariableSizeMeshContainer temp_grid(*this);
    temp_grid.add(extra);
    return temp_grid;
}

template <typename T>
VariableSizeMeshContainer<T>& VariableSizeMeshContainer<T>::operator+=(const VariableSizeMeshContainer<T>& extra)
{
    (*this).add(extra);
    return *this;
}

template <typename T>
T* VariableSizeMeshContainer<T>::operator[](int i)
{
    if (i >= blockNumber)
    {
        std::cerr << "Error: Index out of bounds" << std::endl;
        return nullptr;
    }
    else
    {
        return &V[0] + IA[i];
    }
}

template <typename T>
const T* VariableSizeMeshContainer<T>::operator[](int i) const
{
    if (i >= blockNumber)
    {
        std::cerr << "Error: Index out of bounds" << std::endl;
        return nullptr;
    }
    else
    {
        return &V[0] + IA[i];
    }
}

template <typename T>
int VariableSizeMeshContainer<T>::getBlockNumber() const
{
    return blockNumber;
}

template <typename T>
int VariableSizeMeshContainer<T>::getBlockSize(int i) const
{
    return IA[i+1] - IA[i];
}