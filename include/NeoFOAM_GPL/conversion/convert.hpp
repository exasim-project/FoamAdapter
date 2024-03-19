// SPDX-License-Identifier: GPLv-3.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

NeoFOAM::Vector convert(const Foam::vector& Type)
{
    return NeoFOAM::Vector(Type[0], Type[1], Type[2]);
};

Foam::vector convert(const NeoFOAM::Vector& Type)
{
    return Foam::vector(Type(0), Type(1),Type(2));
};