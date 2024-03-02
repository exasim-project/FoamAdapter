// SPDX-License-Identifier: GPLv-3.0
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

NeoFOAM::vector convert(const Foam::vector& Type)
{
    return NeoFOAM::vector(Type[0], Type[1], Type[2]);
};

Foam::vector convert(const NeoFOAM::vector& Type)
{
    return Foam::vector(Type(0), Type(1),Type(2));
};