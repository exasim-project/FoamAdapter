// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 NeoFOAM authors
#pragma once
#include "dictionary.H"
#include "NeoN/core/dictionary.hpp"

namespace Foam
{


template<typename T>
bool convertEntry(NeoN::Dictionary& neoDict, const Foam::entry& entry)
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
    FatalError.throwExceptions(true);
    FatalIOError.throwExceptions(true);
    try
    {
        entry.get<T>();
    }
    catch (const Foam::IOerror& ioErr)
    {
        FatalError.throwExceptions(false);
        FatalIOError.throwExceptions(false);
        return false;
    }
    catch (const Foam::error& err)
    {
        FatalError.throwExceptions(false);
        FatalIOError.throwExceptions(false);
        return false;
    }
    FatalError.throwExceptions(false);
    FatalIOError.throwExceptions(false);
    return true;
}

NeoN::Dictionary readFoamDictionary(const Foam::dictionary& dict);

} // namespace Foam
