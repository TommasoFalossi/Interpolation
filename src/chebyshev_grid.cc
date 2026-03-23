#include "Interpolation/chebyshev_grid.hh"

namespace Interpolation
{
namespace Chebyshev
{

    StandardGrid::StandardGrid(size_t p){
        _p=p;
        _Dij.resize(p+1, vector_d(p+1, 0.));
        _Dij[0][0] = (2 * p*p +1) / 6.;
        _Dij[p][p] = - _Dij[0][0];


        for(size_t i=0; i<p+1; i++){


            double sign = i % 2 == 0 ? +1 : -1;
            if(i==0 || i==p){
                _betaj.push_back(sign*0.5);
            }
            else{
                _betaj.push_back(sign*1.);
            }
            _tj.push_back(cos(i*M_PI/ static_cast<double>(p)));


            if(i!=0 && i!=p){

                _Dij[i][i] = -0.5 * _tj[i] / (1. -pow(_tj[i], 2));
            }
            
            
        }

        for(size_t i=0; i<p+1; i++){
            for(size_t j =0; j<p+1;j++){
                if(i==j){
                    continue;
                }
                _Dij[i][j] = - (_betaj[i]/_betaj[j]) / (_tj[i]- _tj[j]);
            }

        }
        

    };
} // namespace Chebyshev
} // namespace Interpolation