#include "Interpolation/generic_grid.hh"
#include<stdexcept>

namespace Interpolation
{
namespace Generic
{

    StandardGrid::StandardGrid(const vector_d &input) : _tj(input)
    {

        if(input.size() <= 1){
            throw std::invalid_argument(
                "[Generic::StandardGrid]: input vector cannot be of less than two elements.");

        }

        std::sort(_tj.begin(), _tj.end());

        if(_tj.front() < -1.){
            throw std::invalid_argument("Invalid vector");
        }
        if(_tj.back() > +1.){
            throw std::invalid_argument("Invalid vector");
        }

        if (std::abs(_tj.front() - (-1.)) > 1.0e-12){
            _tj.insert(_tj.begin(), -1.);
        } else {
            _tj.front() = -1.;
        }

        if (std::abs(_tj.back() - 1.) > 1.0e-12){
            _tj.push_back( +1.);
        } else {
            _tj.back() = +1.;
        }

        _p = _tj.size() - 1;
        _lambdaj.resize(_p + 1, 1.);

        for(size_t j = 0; j<=_p; j++){
            for(size_t i=0; i<j; i++){
                _lambdaj[j] *= _tj[j] - _tj[i];
            }
            for(size_t i = j+1; i<= _p; i++){
                _lambdaj[j] *= _tj[j] - _tj[i];
            }

            _lambdaj[j] = 1. / _lambdaj[j];
        }

        //Derivative Matrix
        _Dij.resize(_p + 1, std::vector<double>(_p + 1, 0.));

        for (size_t i = 0; i<= _p; i++){
            //Diagonal Elements
            for(size_t n=0; n<1; n++ ){
                _Dij[i][i] += 1. /(_tj[i] - _tj[n]);
            }
            for(size_t n= i+1; n<=_p; n++){
                _Dij[i][i] += 1. /(_tj[i] - _tj[n]);
            }

            //_Dij[i][j] = 1./ (_tj[i]- _tj[j]);

        }




    };

} // namespace Generic
} // namespace Interpolation