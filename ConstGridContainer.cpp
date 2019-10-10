#include <vector>
#include <iostream>


//Realizing one-dimensional effective-accessible and storing block vector with matrix interface.
template <typename T>
class ConstGridContainer
{

    int M;//blockSize - offset
    std::vector<T> V;//main grid vector

public:

    ConstGridContainer(int blockSize)//number is size of each vector's block
    {
        M = blockSize;
    }


    ConstGridContainer(ConstGridContainer<T>& source)//number is size of each vector's block
    {
        M = source.getBlockSize();
        V.resize(source.getBlockSize()*source.getBlocksNumber());//allocating

        for(int i = 0;i<source.getBlocksNumber();++i)
        {
            for(int j = 0;j<getBlockSize();++j)
            {
                operator[](i)[j] = source[i][j];
            }
        }
    }


    T* operator[](int i)
    {
        return &V[0] + getBlockSize()*i;
    }

    //Can only add elements of vectors which sizes are divided without remainder into {M}
    bool add(const std::vector<T>& extra)
    {
        if (extra.size() % M != 0)
        {
            std::cerr << "Vector size must be a multiple of " <<  M << std::endl;
            return false;
        } else
            {
                V.reserve(extra.size());
                for(typename std::vector<T>::const_iterator it = extra.begin(); it != extra.end();++it)  {
                    V.push_back(*it);
                }
                return true;
            }
    }

    //Add another GridContainer elements to our container
    bool add(ConstGridContainer<T>& extra)
    {
        if (extra.getBlockSize() != getBlockSize())
        {
            std::cerr << "Incompatible block sizes" << std::endl;
            return false;
        } else
            {
                int old_size = getBlocksNumber();
                V.resize((getBlocksNumber() + extra.getBlocksNumber()) * getBlockSize());
                for(int i = old_size;i<getBlocksNumber();++i)
                {
                    for(int j = 0;j<getBlockSize();++j)
                    {
                        operator[](i)[j] = extra[i - old_size][j];
                    }
                }
                return true;
            }
    }


    int getBlocksNumber() const
    {
        return V.size()/getBlockSize();
    }

    int getBlockSize() const
    {
        return M;
    }

    void printContainer(std::ostream& out = std::cout)
    {
        for(int i = 0;i<getBlocksNumber();++i)
        {
            for(int j = 0;j < M;++j)
            {
                out << operator[](i)[j];
                out << " ";
            }
            out << std::endl;
        }
    }


};
