// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once

#include "dictionary.H"
#include "NeoFOAM/core/dictionary.hpp"

namespace Foam
{

// std::vector<std::function<bool(NeoFOAM::Dictionary&, const Foam::entry&)>> mapEntries;

template<typename T>
bool convertEntry(NeoFOAM::Dictionary& neoDict, const Foam::entry& entry)
{
    if (checkEntryType<T>(entry))
    {
        std::string keyword = entry.keyword();
        neoDict.insert(keyword, convert(entry.get<T>()));
        return true;
    }
    return false;
}


template<typename T>
bool checkEntryType(const Foam::entry& entry)
{
    const bool throwingError = FatalError.throwExceptions();
    const bool throwingIOerr = FatalIOError.throwExceptions();
    try
    {
        entry.get<T>();
    }
    catch (const Foam::IOerror& ioErr)
    {
        return false;
    }
    catch (const Foam::error& err)
    {
        return false;
    }

    return true;
}

NeoFOAM::Dictionary readFoamDictionary(const Foam::dictionary& dict);

} // namespace Foam
