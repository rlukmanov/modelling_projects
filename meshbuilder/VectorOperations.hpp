#pragma once
#include <vector>
#include <omp.h>

namespace vops
{
    template<typename VT, typename ST>
    std::vector<VT> axpby(const std::vector<VT>& x,const std::vector<VT>& y,const ST& a,const ST& b,unsigned short threadsNumber = 4)
    {
        std::vector<VT> result;
        size_t size = x.size();

        if (size != y.size())
        {
            std::cerr << "Incompatible sizes!" << std::endl;
            return result;
        }
        else
        {
            result.reserve(size);
            omp_set_num_threads(threadsNumber);
            #pragma omp parallel for 
            for(size_t i = 0; i < size;++i)
            {
                result[i] = a * x[i] + b * y[i];
            }
            return result;
        }
    }

    template<typename VT>
    VT dot(const std::vector<VT>& x,const std::vector<VT>& y,unsigned short threadsNumber = 4)
    {
        VT result = 0;
        size_t size = x.size();

        if (size != y.size())
        {
            std::cerr << "Incompatible sizes!" << std::endl;
            return result;//incorrent
        }
        else
        {
            omp_set_num_threads(threadsNumber);
            #pragma omp parallel for reduction(+:result)
            for(size_t i = 0; i < size;++i)
            {
                result += x[i] * y[i];
            }
            
            return result;
        }
    }
}