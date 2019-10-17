#include <vector>
#include <iostream>
#ifndef __GRID
#define __GRID

//Realizing one-dimensional effective-accessible and storing block vector with matrix interface.
template <typename T>
class GridContainer
{
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

    void inline  printContainer(std::ostream& out = std::cout)
    {
        for(int i = 0;i<getBlocksNumber();++i)
        {
            for(int j = 0;j < getBlockSize(i);++j)
            {
                out << operator[](i)[j];
                out << " ";
            }
            out << std::endl;
        }
    }
};

template <typename T>
int GridContainer<T>::countBlockNumber(int vectorSize, const std::vector<int>& blockSizes)
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
bool GridContainer<T>::checkSizes(const std::vector<T>& source, const std::vector<int>& blockSizes)
{
    int blSizeElems = 0;
    for(typename std::vector<int>::const_iterator it = blockSizes.begin(); it != blockSizes.end();++it)
    {
        blSizeElems += *it;
    }
    return (source.size() == blSizeElems);
}

//default constructor
template <typename T>
GridContainer<T>::GridContainer()
{
    blockNumber = 0;
}


//Creates GridContainer from source vector.
template <typename T>
GridContainer<T>::GridContainer(const std::vector<T>& source, const std::vector<int>& blockSizes)//number is size of each vector's block
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
GridContainer<T>::GridContainer(GridContainer<T>& source)//number is size of each vector's block
{
    blockNumber = source.getBlocksNumber();
    V.resize(source.getIa()[source.getIa().size() - 1]);//allocating
    IA.resize(source.getIa().size());

    for(int i = 0;i<source.getIa().size();++i)
    {
        IA[i] = source.getIa()[i];
    }

    for(int i = 0;i<source.getBlocksNumber();++i)
    {
        for(int j = 0;j<source.getBlockSize(i);++j)
        {
            operator[](i)[j] = source[i][j];
        }
    }

}

//Can only add elements of vectors which sizes are divided without remainder into {M}
template <typename T>
bool GridContainer<T>::add(const std::vector<T>& extra,const std::vector<int>& blockSizes)
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


//Add another GridContainer elements to our container
template <typename T>
bool GridContainer<T>::add(GridContainer<T>& extra)
{
    int oldBlocksNum = blockNumber;
    int oldOffsetNum = IA.size();
    int offset = IA[IA.size() - 1];

    this->blockNumber += extra.getBlocksNumber();
    IA.resize(IA.size() + extra.getIa().size() - 1);
    V.resize(V.size() + extra.getIa()[extra.getIa().size() - 1]);


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
GridContainer<T> GridContainer<T>::operator+(GridContainer<T>& extra)
{
    GridContainer temp_grid(*this);
    temp_grid.add(extra);
    return temp_grid;
}

template <typename T>
GridContainer<T>& GridContainer<T>::operator+=(GridContainer<T>& extra)
{
    (*this).add(extra);
    return *this;
}

template <typename T>
T* GridContainer<T>::operator[](int i)
{
    if (i >= getBlocksNumber())
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
int GridContainer<T>::getBlocksNumber() const
{
    return blockNumber;
}

template <typename T>
int GridContainer<T>::getBlockSize(int i) const
{
    return IA[i+1] - IA[i];
}

template <typename T>
const std::vector<int>&  GridContainer<T>::getIa() const
{
    return IA;
}

#endif