
from framework import model, step, FixedOuter, Simulation, attach_refs_and_log, run_time_step

@model("DomainA")
class DomainA:
    controller = FixedOuter(iterations=1)
    @step(order=10)
    def interact(self, ctrl):
        # Access DomainB via refs
        b = self.refs["b"]
        self.log.append(("A", "interact", b.__class__.__name__))

@model("DomainB")
class DomainB:
    controller = FixedOuter(iterations=1)
    @step(order=10)
    def respond(self, ctrl):
        self.log.append(("B", "respond"))

def test_cross_domain_coupling():
    config = {
        "regions": [{"name": "a"}, {"name": "b"}],
        "solvers": {"a": "DomainA", "b": "DomainB"},
    }
    sim = Simulation(config)
    attach_refs_and_log(sim)
    run_time_step(sim)
    a_log = sim.domains["a"].model.log
    b_log = sim.domains["b"].model.log
    # DomainA should log interaction with DomainB
    assert a_log == [("A", "interact", "DomainB")]
    assert b_log == [("B", "respond")]
