// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2023 FoamAdapter authors
#pragma once

#include "NeoN/NeoN.hpp"

#include "fvMesh.H"
#include "volFields.H"

#include "FoamAdapter/auxiliary/convert.hpp"

namespace FoamAdapter
{

namespace detail
{

template<class DestField, class SrcField>
void copy_impl(DestField& dest, const SrcField src)
{
    NF_ASSERT_EQUAL(dest.size(), src.size());
    auto srcHost = src.copyToHost();
    auto srcView = srcHost.view();
    for (int i = 0; i < dest.size(); i++)
    {
        dest[i] = convert(srcView[i]);
    }
}
}

void write(NeoN::scalarVector& sf, const Foam::fvMesh& mesh, const std::string fieldName);

void write(NeoN::Vector<NeoN::Vec3>& sf, const Foam::fvMesh& mesh, const std::string fieldName);

} // namespace Foam
