import framework
import my_model

CONFIG = {
    "regions": [{"name": "fluid"}, {"name": "solid"}],
    "solvers": {"fluid": "singlePhasePIMPLE", "solid": "solidConduction"},
}

def test_pressure_repeats_and_momentum_once():
    sim = framework.Simulation(CONFIG)
    framework.attach_refs_and_log(sim)
    framework.run_time_step(sim)
    fluid_log = sim.domains["fluid"].model.log
    # Pressure repeats 3x per outer (2 outer)
    pressure = [x for x in fluid_log if x[1] == "pressure"]
    assert len(pressure) == 6
    # Momentum runs once per outer (2 outer)
    momentum = [x for x in fluid_log if x[1] == "momentum"]
    assert len(momentum) == 2
