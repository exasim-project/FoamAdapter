// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors

#include "FoamAdapter/conversion/convert.hpp"

namespace Foam
{

NeoFOAM::Vector convert(const Foam::vector& in) { return NeoFOAM::Vector(in[0], in[1], in[2]); };

NeoFOAM::scalar convert(const Foam::scalar& in) { return in; };

NeoFOAM::label convert(const Foam::label& in) { return in; };

std::string convert(const Foam::word& type) { return type; };

NeoFOAM::TokenList convert(const Foam::ITstream& stream)
{
    NeoFOAM::TokenList tokens {};
    for (const auto& token : stream)
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
Foam::vector convert(const NeoFOAM::Vector& in) { return Foam::vector(in(0), in(1), in(2)); };

}; // namespace Foam
