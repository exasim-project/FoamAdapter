import my_model
from framework import list_models

def test_registry_contains_models():
    models = list_models()
    assert "singlePhasePIMPLE" in models
    assert "solidConduction" in models


from framework import model, step, FixedOuter, Simulation, attach_refs_and_log, run_time_step

@model("customSolver")
class CustomSolver:
    controller = FixedOuter(iterations=1)
    @step(order=10)
    def run(self, ctrl):
        self.log.append(("custom", ctrl.outer_iter, ctrl.correct_iter))

def test_custom_solver_registration_and_run():
    config = {
        "regions": [{"name": "custom"}],
        "solvers": {"custom": "customSolver"},
    }
    sim = Simulation(config)
    attach_refs_and_log(sim)
    run_time_step(sim)
    log = sim.domains["custom"].model.log
    assert log == [("custom", 0, 0)]
