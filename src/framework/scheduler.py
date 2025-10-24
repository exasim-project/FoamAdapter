def collect_steps(model):
    # Use StepMeta dataclasses from model._steps
    steps = getattr(model.__class__, "_steps", [])
    # Sort by order
    steps = sorted(steps, key=lambda s: s.order)
    # Return as (order, repeat, func) tuples for compatibility
    return [(s.order, s.repeat, s.func) for s in steps]

def attach_refs_and_log(sim):
    # Attach refs and log to each model
    for name, domain in sim.domains.items():
        domain.model.refs = {k: v.model for k, v in sim.domains.items() if k != name}
        domain.model.log = []

def run_time_step(sim):
    for name, domain in sim.domains.items():
        model = domain.model
        controller = getattr(model, "controller")
        steps = collect_steps(model)
        for ctrl in controller():
            for order, repeat, fn in steps:
                bound_fn = getattr(model, fn.__name__)
                if repeat == "correct":
                    bound_fn(ctrl)
                else:
                    if ctrl.correct_iter == 0:
                        bound_fn(ctrl)
