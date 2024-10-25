// SPDX-License-Identifier: GPL-3.0-or-later
// SPDX-FileCopyrightText: 2024 NeoFOAM authors

NeoFOAM::Executor exec = createExecutor(runTime.controlDict());
std::unique_ptr<MeshAdapter> meshPtr = createMesh(exec, runTime);
MeshAdapter& mesh = *meshPtr;
