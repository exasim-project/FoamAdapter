// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 NeoFOAM authors

#include "FoamAdapter/conversion/convert.hpp"
#include "FoamAdapter/readers/foamDictionary.hpp"
#include "vector.H"
#include <functional>
namespace Foam
{


static std::vector<std::function<bool(NeoFOAM::Dictionary&, const Foam::entry&)>> mapEntries =
    {[](NeoFOAM::Dictionary& neoDict, const Foam::entry& entry)
     {
         if (checkEntryType<Foam::label>(entry))
         {
             Foam::token::tokenType tt = entry.stream()[0].type();
             if (tt == Foam::token::tokenType::LABEL)
             {
                 neoDict.insert(entry.keyword(), convert(entry.get<Foam::label>()));
             }
             else
             {
                 neoDict.insert(entry.keyword(), entry.get<Foam::scalar>());
             }
             return true;
         }
         return false;
     },
     &convertEntry<Foam::scalar>,
     &convertEntry<Foam::vector>,
     &convertEntry<Foam::word>};


void insertEntry(NeoFOAM::Dictionary& neoDict, const Foam::entry& entry)
{
    std::string keyword = entry.keyword();
    for (auto& mapEntry : mapEntries)
    {
        if (mapEntry(neoDict, entry))
        {
            return;
        }
    }
    std::string entryString;
    if (entry.isStream())
    {
        entryString = entry.stream().toString();
    }
    throw std::runtime_error(
        "No known conversion for the key: " + keyword
        + " \n"
          "and the following entry: "
        + entryString
    );
}

void readFoamDictionary(const Foam::dictionary& dict, NeoFOAM::Dictionary& neoDict)
{
    for (auto& entry : dict)
    {
        std::string keyword = entry.keyword();
        if (entry.isDict())
        {
            NeoFOAM::Dictionary subDict;
            readFoamDictionary(entry.dict(), subDict);
            neoDict.insert(entry.keyword(), subDict);
        }
        else
        {
            insertEntry(neoDict, entry);
        }
    }
}

NeoFOAM::Dictionary readFoamDictionary(const Foam::dictionary& dict)
{
    NeoFOAM::Dictionary neoDict;
    readFoamDictionary(dict, neoDict);
    return neoDict;
}

} // namespace Foam
