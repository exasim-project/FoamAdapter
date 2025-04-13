// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "volFields.H"

#include "NeoN/core/primitives/label.hpp"
#include "NeoN/core//primitives/scalar.hpp"
#include "NeoN/core/primitives/vector.hpp"
#include "NeoN/core/tokenList.hpp"
#include "volFields.H"
// #include "ITstream.H"

namespace Foam
{

NeoN::Vector convert(const Foam::vector& Type);
NeoN::scalar convert(const Foam::scalar& Type);
std::string convert(const Foam::word& Type);

NeoN::TokenList convert(const Foam::ITstream& Type);

NeoN::label convert(const Foam::label& Type);

// To Foam
Foam::vector convert(const NeoN::Vector& Type);

}; // namespace Foam
