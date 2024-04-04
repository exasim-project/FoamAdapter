// SPDX-License-Identifier: GPLv-3.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#include "FoamAdapter/conversion/convert.hpp"


namespace Foam
{

NeoFOAM::Vector convert(const Foam::vector& Type)
{
    return NeoFOAM::Vector(Type[0], Type[1], Type[2]);
};

NeoFOAM::scalar convert(const Foam::scalar& Type)
{
    return Type;
};

NeoFOAM::label convert(const Foam::label& Type)
{
    return Type;
};

// To Foam
Foam::vector convert(const NeoFOAM::Vector& Type)
{
    return Foam::vector(Type(0), Type(1),Type(2));
};

}; // namespace Foam
