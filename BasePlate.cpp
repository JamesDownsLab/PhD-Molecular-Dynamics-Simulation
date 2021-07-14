//
// Created by ppxjd3 on 14/07/2021.
//

#include "BasePlate.h"

void BasePlate::update(double Time) {
    _z = _A * sin(_omega * Time);
    _vz = _A * _omega * cos(_omega * Time);
}
