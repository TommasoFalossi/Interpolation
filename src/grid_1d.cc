#include "Interpolation/grid_1d.hh"
#include <algorithm>

namespace Interpolation{


    Grid1D::Grid1D(const SingleDiscretizationInfo &d_info){

        _d_info = d_info;

        for(size_t i= 0; i < _d_info.intervals.size(); i++){
            _stored_grids.insert({_d_info.grid_sizes[i], Chebyshev::StandardGrid(_d_info.grid_sizes[i])});
        }
        
        ////////////////////////////////////////////////////////

        for(size_t i = 0; i<_d_info.intervals.size(); i++){

            vector_d t_j = _stored_grids[_d_info.grid_sizes[i]]._tj;
            std::reverse(t_j.begin(), t_j.end());
            double a_inter = _d_info.intervals[i].first;
            double b_inter = _d_info.intervals[i].second;

            for(double x : t_j){
                _coord_inter.push_back(a_inter + (b_inter - a_inter)*(x+1.0)/2.0);
            }

        }

        for(size_t i=0; i<_coord_inter.size(); i++){
            _coord.push_back(_d_info.to_phys_space(_coord_inter[i]));
        }


        ////////////////////////////////////////////////////////

    }


}

