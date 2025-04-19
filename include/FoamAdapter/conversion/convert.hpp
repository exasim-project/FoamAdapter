// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "volFields.H"

#include "NeoN/NeoN.hpp"
// #include "ITstream.H"

namespace Foam
{

NeoN::Vec3 convert(const Foam::vector& Type);

NeoN::scalar convert(const Foam::scalar& Type);

std::string convert(const Foam::word& Type);

NeoN::TokenList convert(const Foam::ITstream& Type);

NeoN::label convert(const Foam::label& Type);

// To Foam
Foam::vector convert(const NeoN::Vec3& Type);

}; // namespace Foam
