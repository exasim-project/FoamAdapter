import typer
import sys

app = typer.Typer()


@app.callback()
def cli():
    """FoamAdapter CLI"""
    pass


# Solver command group
solver_app = typer.Typer()
app.add_typer(solver_app, name="solver")


@solver_app.command(
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True}
)
def icofoam(
    ctx: typer.Context,
    check_inputs: bool = typer.Option(
        False, "--check_inputs", help="Check inputs and mesh before running"
    ),
):
    """Run the icoFoam solver."""
    from foamadapter.solver.icoFoam import main as icofoam_main

    # Only pass the extra args (not the Typer command path)
    argv = [sys.argv[0]] + [str(arg) for arg in ctx.args]
    icofoam_main(argv)


@solver_app.command(
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True}
)
def pimplefoam(
    ctx: typer.Context,
    check_inputs: bool = typer.Option(
        False, "--check_inputs", help="Check inputs and mesh before running"
    ),
):
    """Run the pimpleFoam solver."""

    from foamadapter.solver.pimpleFoam import main as pimplefoam_main

    # Only pass the extra args (not the Typer command path)
    if check_inputs:
        # Here you would implement the logic to check inputs and mesh
        print("Checking inputs and mesh... (not implemented in this example)")
        print("You can implement this check in the solver code.", ctx.args, "asd ", sys.argv)
        return
    argv = [sys.argv[0]] + [str(arg) for arg in ctx.args]
    pimplefoam_main(argv)
