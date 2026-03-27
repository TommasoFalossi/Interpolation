#include "Interpolation/chebyshev_grid.hh"
#include<stdexcept>
#include<iostream>

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

    double StandardGrid::poli_weight(double t, size_t j) const{

        if(std::abs(t - this->_tj[j]) < 1.0e-15) {
                return 1.;
            }

        double den = 0;
        for(size_t i=0; i< this->_betaj.size(); i++){
            if(std::abs(t - this->_tj[i]) < 1.0e-15) {
                return 0.;
            }
            den += (this->_betaj[i]/(t-this->_tj[i]));
        }

        double b_i=this->_betaj[j]/(t-this->_tj[j])/den;

        return b_i;
    };

    double StandardGrid::poli_weight(double t, size_t j, double den) const{

        if(std::abs(t - this->_tj[j]) < 1.0e-15) {
                return 1.;
            }

        double b_i=this->_betaj[j]/(t-this->_tj[j])/den;

        return b_i;

    };

    double StandardGrid::interpolate(double t, const vector_d &fj, size_t start, size_t end) const{

        if(t<-1 || t>1){
            throw std::domain_error("StandardGrid::interpolate: t must be in [-1, 1]");
        }
        if(end-start != _p){
            throw std::domain_error("StandardGrid::interpolate: end-start must be equal to p");
        }

        /*
        double p = 0;
        for(size_t i=0; i<=this->_p; i++){

            double b=this->poli_weight(t, i);

            p+=fj[i+start]*b;
        }

        return p;
        */

        double den = 0;
        for(size_t i=0; i< this->_betaj.size(); i++){
            if(std::abs(t - this->_tj[i]) < 1.0e-15) {
                return fj[i+start];
            }
            den += (this->_betaj[i]/(t-this->_tj[i]));
        }

        double p = 0;
        for(size_t i=0; i<=this->_p; i++){

            p += fj[i+start] * this->poli_weight(t, i, den);
        }

        return p;

    };

    vector_d StandardGrid::discretize(const std::function<double(double)> &fnc) const{

        vector_d fj(_p + 1, 0.);
        for(size_t i=0; i<=_p; i++){
            fj[i]=fnc(_tj[i]);
        }

        return fj;
    };


    double StandardGrid::interpolate_der(double t, const vector_d &fj, size_t start, size_t end) const{

        if(t<-1 || t>1){
            throw std::domain_error("StandardGrid::interpolate: t must be in [-1, 1]");
        }
        if(end-start != _p){
            throw std::domain_error("StandardGrid::interpolate: end-start must be equal to p");
        }


        double p = 0;
        for(size_t i=start, j=0; i<=end; i++, j++){
            p+=fj[i]*poli_weight_der(t, j);
        }

        return p;
    };

    double StandardGrid::poli_weight_der(double t, size_t j) const{
        double l=0;
        for(size_t i = 0 ; i <= _p; i++){
            l+= poli_weight(t, i)* _Dij[j][i];
        }

        return l;
    };

    double StandardGrid::poli_weight_der(double t, size_t j, double den) const{
        double l=0;
        for(size_t i = 0 ; i <= _p; i++){
            l+= poli_weight(t, i, den)* _Dij[j][i];
        }

        return l;
    };

    void StandardGrid::apply_D(vector_d &fj, size_t start, size_t end) const{

        vector_d f_tilde = fj;

        for(size_t i=start; i<=end; i++){
            for(size_t j=start; j<=end; j++){
                fj[i] = 0;
            }
        }


    }

} // namespace Chebyshev
} // namespace Interpolation