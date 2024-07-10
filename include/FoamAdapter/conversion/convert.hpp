// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "NeoFOAM/primitives/label.hpp"
#include "NeoFOAM/primitives/scalar.hpp"
#include "NeoFOAM/primitives/vector.hpp"
#include "volFields.H"

namespace Foam
{

NeoFOAM::Vector convert(const Foam::vector& Type);
NeoFOAM::scalar convert(const Foam::scalar& Type);

NeoFOAM::label convert(const Foam::label& Type);

// To Foam
Foam::vector convert(const NeoFOAM::Vector& Type);

}; // namespace Foam
