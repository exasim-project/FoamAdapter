// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2023 FoamAdapter authors

#pragma once

#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <Kokkos_Core.hpp>

#include "NeoN/NeoN.hpp"

#include <string>
#include <utility>

// This class shows how to implement a simple generator for Catch tests
class ExecutorGenerator final :
    public Catch::Generators::IGenerator<std::pair<std::string, NeoN::Executor>>
{
public:

    std::size_t i = 0;
    std::vector<std::pair<std::string, NeoN::Executor>> execs {};
    std::pair<std::string, NeoN::Executor> current_exec = {
        "SerialExecutor",
        NeoN::SerialExecutor {}
    };

    ExecutorGenerator()
    {
#if defined(KOKKOS_ENABLE_OPENMP)
        execs.push_back({"ParallelCPUExecutor", NeoN::CPUExecutor {}});
#elif defined(KOKKOS_ENABLE_THREADS)
        execs.push_back({"ParallelCPUExecutor", NeoN::CPUExecutor {}});
#endif
#if defined(KOKKOS_ENABLE_CUDA)
        execs.push_back({"GPUExecutor", NeoN::GPUExecutor {}});
#elif defined(KOKKOS_ENABLE_HIP)
        execs.push_back({"GPUExecutor", NeoN::GPUExecutor {}});
#elif defined(KOKKOS_ENABLE_SYCL)
        execs.push_back({"GPUExecutor", NeoN::GPUExecutor {}});
#endif
    }

    std::pair<std::string, NeoN::Executor> const& get() const override;
    bool next() override
    {
        if (i >= execs.size()) return false;
        current_exec = execs[i];
        i++;
        return true;
    }
};

// Avoids -Wweak-vtables
std::pair<std::string, NeoN::Executor> const& ExecutorGenerator::get() const
{
    return current_exec;
}

Catch::Generators::GeneratorWrapper<std::pair<std::string, NeoN::Executor>> allAvailableExecutor()
{
    return Catch::Generators::GeneratorWrapper<std::pair<std::string, NeoN::Executor>>(
        Catch::Detail::make_unique<ExecutorGenerator>()
    );
}
