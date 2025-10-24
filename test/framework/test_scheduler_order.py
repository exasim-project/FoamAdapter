
from framework import model, step, FixedOuter, Simulation, attach_refs_and_log, run_time_step

@model("TestFluid")
class TestFluid:
    controller = FixedOuter(iterations=2)
    @step(order=10)
    def momentum(self, ctrl):
        self.log.append(("fluid", "momentum", ctrl.outer_iter, ctrl.correct_iter))
    @step(order=20)
    def energy(self, ctrl):
        self.log.append(("fluid", "energy", ctrl.outer_iter, ctrl.correct_iter))

@model("TestSolid")
class TestSolid:
    controller = FixedOuter(iterations=1)
    @step(order=10)
    def conduction(self, ctrl):
        self.log.append(("solid", "conduction", ctrl.outer_iter, ctrl.correct_iter))

def test_fluid_and_solid_step_order():
    config = {
        "regions": [{"name": "fluid"}, {"name": "solid"}],
        "solvers": {"fluid": "TestFluid", "solid": "TestSolid"},
    }
    sim = Simulation(config)
    attach_refs_and_log(sim)
    run_time_step(sim)
    fluid_log = sim.domains["fluid"].model.log
    solid_log = sim.domains["solid"].model.log
    # Fluid: momentum before energy, for each outer iteration
    fluid_steps = [x[1] for x in fluid_log if x[1] in ("momentum", "energy")]
    assert fluid_steps == ["momentum", "energy", "momentum", "energy"]
    # Solid ran once
    assert solid_log == [("solid", "conduction", 0, 0)]
