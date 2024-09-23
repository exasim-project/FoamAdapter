// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "volFields.H"

#include "NeoFOAM/core/primitives/label.hpp"
#include "NeoFOAM/core//primitives/scalar.hpp"
#include "NeoFOAM/core/primitives/vector.hpp"

namespace Foam
{

NeoFOAM::Vector convert(const Foam::vector& Type);
NeoFOAM::scalar convert(const Foam::scalar& Type);

NeoFOAM::label convert(const Foam::label& Type);

// To Foam
Foam::vector convert(const NeoFOAM::Vector& Type);

}; // namespace Foam
