// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "FoamAdapter/conversion/convert.hpp"

namespace Foam
{

NeoFOAM::Vector convert(const Foam::vector& in) { return NeoFOAM::Vector(in[0], in[1], in[2]); };

NeoFOAM::scalar convert(const Foam::scalar& in) { return in; };

NeoFOAM::label convert(const Foam::label& in) { return in; };

// To Foam
Foam::vector convert(const NeoFOAM::Vector& in) { return Foam::vector(in(0), in(1), in(2)); };

}; // namespace Foam
