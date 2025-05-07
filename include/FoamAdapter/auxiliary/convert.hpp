// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 FoamAdapter authors
#pragma once

#include "NeoN/NeoN.hpp"

#include "fvc.H"

namespace FoamAdapter
{
  // To NeoN
NeoN::Vec3 convert(const Foam::vector& Type);

NeoN::scalar convert(const Foam::scalar& Type);

std::string convert(const Foam::word& Type);

NeoN::TokenList convert(const Foam::ITstream& Type);

NeoN::label convert(const Foam::label& Type);

NeoN::Dictionary convert(const Foam::dictionary& dict);

template<typename T>
bool checkEntryType(const Foam::entry& entry)
{
    Foam::FatalError.throwExceptions(true);
    Foam::FatalIOError.throwExceptions(true);
    try
    {
        entry.get<T>();
    }
    catch (const Foam::IOerror& ioErr)
    {
        Foam::FatalError.throwExceptions(false);
        Foam::FatalIOError.throwExceptions(false);
        return false;
    }
    catch (const Foam::error& err)
    {
        Foam::FatalError.throwExceptions(false);
        Foam::FatalIOError.throwExceptions(false);
        return false;
    }
    Foam::FatalError.throwExceptions(false);
    Foam::FatalIOError.throwExceptions(false);
    return true;
}

template<typename T>
bool insert(NeoN::Dictionary& neoDict, const Foam::entry& entry)
{
    if (checkEntryType<T>(entry))
    {
        std::string keyword = entry.keyword();
        neoDict.insert(keyword, convert(entry.get<T>()));
        return true;
    }
    return false;
}

NeoN::TokenList convert(const Foam::ITstream& stream);

// To Foam
Foam::vector convert(const NeoN::Vec3& Type);

}; // namespace FoamAdapter
