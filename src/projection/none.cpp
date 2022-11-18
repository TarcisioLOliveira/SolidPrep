/*
 *   Copyright (C) 2022 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   SolidPrep is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   SolidPrep is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with SolidPrep.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#include <cblas.h>
#include "projection/none.hpp"

namespace projection{

void None::update(const size_t iteration){
    (void)iteration;
}

void None::project_densities(const std::vector<double>& x, std::vector<double>& new_x) const{
    cblas_dcopy(x.size(), x.data(), 1, new_x.data(), 1);
}

void None::project_gradient(const std::vector<double>& df, std::vector<double>& new_df) const{
    cblas_dcopy(df.size(), df.data(), 1, new_df.data(), 1);
}

}
