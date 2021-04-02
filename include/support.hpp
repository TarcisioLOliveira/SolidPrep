/*
 *   Copyright (C) 2021 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of SolidPrep
 *
 *   Foobar is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Foobar is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with SolidPrep.  If not, see <https://www.gnu.org/licenses/>.
 *
 */
#ifndef SUPPORT_HPP
#define SUPPORT_HPP

#include <memory>
#include <vector>
#include "cross_section.hpp"
#include "BRepClass3d_SolidClassifier.hxx"

class Support{
    public:

    /**
     * Creates a Support object for 2D problems.
     *
     * @param X Whether it resists displacement along X.
     * @param Y Whether it resists displacement along Y.
     * @param MZ Whether it resists bending perpendicular to Z.
     * @param cross_section Maximum dimensions of the supporting face.
     */
    Support(bool X, bool Y, bool MZ, CrossSection cross_section);

    /**
     * Creates a Support object for 3D problems.
     *
     * @param X Whether it resists displacement along X.
     * @param Y Whether it resists displacement along Y.
     * @param Z Whether it resists displacement along Z.
     * @param MX Whether it resists bending perpendicular to X.
     * @param MY Whether it resists bending perpendicular to Y.
     * @param MZ Whether it resists bending perpendicular to Z.
     * @param cross_section Maximum dimensions of the supporting face.
     */
    Support(bool X, bool Y, bool Z, bool MX, bool MY, bool MZ, CrossSection cross_section);

    inline int F_DOF() const{
        return this->fdof;
    }
    inline int M_DOF() const{
        return this->mdof;
    }
    inline int DOF() const{
        return this->fdof + this->mdof;
    }

    const bool X, Y, Z, MX, MY, MZ;
    CrossSection S;
    private:
    int fdof, mdof;
};

#endif
