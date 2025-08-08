import typer

app = typer.Typer()

@app.callback()
def cli():
    """FoamAdapter CLI"""
    pass

@app.command()

def hello(name: str = "World"):
    """Say hello."""
    typer.echo(f"Hello {name}!")



# Solver command group
solver_app = typer.Typer()
app.add_typer(solver_app, name="solver")

@solver_app.command(context_settings={"allow_extra_args": True, "ignore_unknown_options": True})
def icofoam(ctx: typer.Context):
    """Run the icoFoam solver."""
    import sys
    from foamadapter.cli.solver.icoFoam import main as icofoam_main
    # Only pass the extra args (not the Typer command path)
    argv = [sys.argv[0]] + [str(arg) for arg in ctx.args]
    icofoam_main(argv)
