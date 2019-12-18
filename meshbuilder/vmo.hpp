#pragma once

#include <vector>
#include <omp.h>
#include "Sparse.hpp"

//vector matrix operations: consider we have same template arguments
namespace vmo
{
    template<typename VT, typename ST>
    std::vector<VT> axpby(const std::vector<VT>& x,const std::vector<VT>& y,const ST& a,const ST& b,unsigned short threadsNumber = 4)
    {
        std::vector<VT> result;
        size_t size = x.size();
        if (size != y.size())
        {
            std::cerr << "Incompatible sizes in axpby!" << std::endl;
            return result;
        }
        else
        {
            result.resize(size);

            #pragma omp parallel for num_threads(threadsNumber) 
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
            std::cerr << "Incompatible sizes in dot!" << std::endl;
            return result;//incorrent
        }
        else
        {
            #pragma omp parallel for num_threads(threadsNumber) reduction(+:result)
            for(size_t i = 0; i < size;++i)
            {
                result += x[i] * y[i];
            }
            return result;
        }
    }

    template<typename VT>
    std::vector<VT> multiply(const std::vector<VT>& x,VT scalar,unsigned short threadsNumber = 4)
    {
        std::vector<VT> result;

        result.resize(x.size());

        #pragma omp parallel for num_threads(threadsNumber) 
        for(size_t i = 0; i < x.size();++i)
        {
            result[i] = scalar * x[i];
        }
        return result;
    }

    //A = A^(T) > 0
    template<typename M,typename V>
    std::vector<double> conGradSolver(Sparse<M>& A, std::vector<V> b,double eps = 0.001)
    {
        std::vector<double> x_prev(A.getDenseColumns());
        if (A.getDenseColumns() != b.size())
        {
            std::cerr << "Incompatible sizes in Solver!" << std::endl;
        }
        else
        {
            size_t size = A.getDenseColumns();  

            double** rev_M = new double*[size];
            for(size_t i = 0; i < size;++i)
            {
                rev_M[i] = new double[size];
            }

            M** dense = A.getDenseMatrix();

            //forming preconditioner Matrix

            for(size_t i = 0;i < size;++i)
            {
                for(size_t j = 0;j < size;++j)
                {
                    rev_M[i][j] = 0.0;
                }
                rev_M[i][i] = 1 / (dense[i][i] + eps);
            }
            
            SparseELL<double> reverse_M(rev_M,size,size);//preconditioner

            for(size_t i = 0;i < size;++i)
            {
                delete[] dense[i];
                delete[] rev_M[i];
            }   
            delete[] dense;
            delete[] rev_M;
            
            //initializang parameters for algorithm
            
            std::vector<double> x_cur;

            std::vector<double> r_cur;
            std::vector<double> r_prev = b;//we think that x0 = (0,0,0,0,0,0 ... 0)^(T) -> b - Ax = b
            
            std::vector<double> p_cur;//p and p+1 are conjugated relative to A
            std::vector<double> p_prev;

            std::vector<double> q;
            std::vector<double> z;

            double delta_cur,delta_prev;
            double alpha,beta;

            size_t k = 1;
            
            do 
            {
                std::cout << delta_cur << std::endl;
                z = reverse_M.spmv(r_prev);
                delta_cur = dot(r_prev,z);
                if (k == 1)
                {
                    p_cur = z;
                }
                else
                {
                    beta = delta_cur / delta_prev;
                    p_cur = axpby(z,multiply(p_prev,beta),1,1);
                }
                q = A.spmv(p_cur);
                alpha = delta_cur / dot(p_cur,q);

                x_cur = axpby(x_prev,multiply(p_cur,alpha),1,1);

                r_cur = axpby(r_prev,multiply(q,alpha),1,-1);
                if (delta_cur < eps)
                {
                    break;
                }
                
                r_prev = r_cur;
                x_prev = x_cur;
                p_prev = p_cur;
                delta_prev = delta_cur;

                k += 1;
                
            } while(true);

        }
        
        return x_prev;
    }
}