//
// Created by ppxjd3 on 14/07/2021.
//

#include "Particle.h"

std::istream &operator>>(std::istream &is, Particle &p) {
    is >> p.rtd0 >> p._r >> p._m >> p._youngs_modulus >> p._poisson >> p._coeff_res >> p._coeff_fric;
    return is;
}
