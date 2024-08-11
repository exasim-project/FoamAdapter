// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#include "FoamAdapter/conversion/convert.hpp"


namespace Foam
{

NeoFOAM::Vector convert(const Foam::vector& Type)
{
    return NeoFOAM::Vector(Type[0], Type[1], Type[2]);
};

NeoFOAM::scalar convert(const Foam::scalar& Type) { return Type; };

NeoFOAM::label convert(const Foam::label& Type) { return Type; };

std::string convert(const Foam::word& Type) { return Type; };

NeoFOAM::TokenList convert(const Foam::ITstream& Type)
{
    NeoFOAM::TokenList tokens;
    for (const auto& token : Type)
    {
        if (token.isBool())
        {
            tokens.insert(token.boolToken());
            continue;
        }
        if (token.isLabel())
        {
            tokens.insert(token.labelToken());
            continue;
        }
        if (token.isScalar())
        {
            tokens.insert(token.scalarToken());
            continue;
        }
        if (token.isWord())
        {
            tokens.insert(std::string(token.wordToken()));
            continue;
        }
        std::runtime_error("Unsupported token type");
    }
    return tokens;
};

// To Foam
Foam::vector convert(const NeoFOAM::Vector& Type)
{
    return Foam::vector(Type(0), Type(1), Type(2));
};

}; // namespace Foam
