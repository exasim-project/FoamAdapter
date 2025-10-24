from framework import FixedOuter, Simulation, attach_refs_and_log, run_time_step

def test_model_connectivity():
    @model("A")
    class A:
        controller = FixedOuter(iterations=1)
        @step(order=1)
        def connect(self, ctrl):
            b = self.refs["b"]
            self.log.append(("A", "connect", b.__class__.__name__))

    @model("B")
    class B:
        controller = FixedOuter(iterations=1)
        @step(order=1)
        def respond(self, ctrl):
            self.log.append(("B", "respond"))

    config = {
        "regions": [{"name": "a"}, {"name": "b"}],
        "solvers": {"a": "A", "b": "B"},
    }
    sim = Simulation(config)
    attach_refs_and_log(sim)
    run_time_step(sim)
    a_log = sim.domains["a"].model.log
    b_log = sim.domains["b"].model.log
    assert a_log == [("A", "connect", "B")]
    assert b_log == [("B", "respond")]
from framework import model, step, get_model

def test_model_and_step_decorators():
    @model("TestModel")
    class TestModel:
        @step(order=1)
        def foo(self, ctrl):
            pass

        @step(order=2, repeat="correct")
        def bar(self, ctrl):
            pass

    # Check model registration
    assert get_model("TestModel") is TestModel

    # Check step metadata via StepMeta dataclass
    foo_meta = TestModel.foo._step_meta
    bar_meta = TestModel.bar._step_meta
    assert foo_meta.name == "foo"
    assert foo_meta.order == 1
    assert foo_meta.repeat is None
    assert bar_meta.name == "bar"
    assert bar_meta.order == 2
    assert bar_meta.repeat == "correct"

    # Check number of steps using class method
    assert TestModel.count_steps() == 2
