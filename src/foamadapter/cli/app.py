import typer
import sys
from pydantic import ValidationError, BaseModel, Field
from typing import List, Optional
from foamadapter.turbulence.incompressible import TurbulenceModel

app = typer.Typer()


@app.callback()
def cli():
    """FoamAdapter CLI"""
    pass


# Dummy Pydantic model for demonstration
class SolverConfig(BaseModel):
    start_time: float = Field(gt=0, description="Start time must be positive")
    end_time: float = Field(gt=0, description="End time must be positive")
    time_step: float = Field(
        gt=0, le=1, description="Time step must be between 0 and 1"
    )
    mesh_quality: Optional[float] = Field(
        None, ge=0, le=1, description="Mesh quality score"
    )
    boundary_conditions: List[str] = Field(
        min_items=1, description="At least one boundary condition required"
    )


def format_pydantic_errors(validation_error: ValidationError) -> str:
    """Format Pydantic validation errors in a user-friendly way."""
    formatted_errors = []

    for error in validation_error.errors():
        field_path = " -> ".join(str(loc) for loc in error["loc"])
        msg = error["msg"]
        formatted_errors.append(
            f"  • {typer.style(field_path, fg=typer.colors.CYAN, bold=True)}: {msg}"
        )

    return "\n".join(formatted_errors)


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
    case_dir: str = typer.Option(".", "--case", help="Case directory path"),
):
    """Run the pimpleFoam solver."""

    from foamadapter.solver.pimpleFoam import PimpleFoam

    # Only pass the extra args (not the Typer command path)
    argv = [sys.argv[0]] + [str(arg) for arg in ctx.args]
    
    if check_inputs:
        try:
            # Get the registry from PimpleFoam
            registry = PimpleFoam(argv).inputs()
            
            # Validate the case using the registry
            is_valid, errors = registry.validate_case(case_dir)
            
            if not is_valid:
                typer.echo(
                    typer.style(
                        "❌ Input validation failed:", fg=typer.colors.RED, bold=True
                    )
                )
                for error in errors:
                    typer.echo(f"  • {error}")
                typer.echo(
                    typer.style(
                        "Fix the above errors before running the solver.",
                        fg=typer.colors.YELLOW,
                    )
                )
                raise typer.Exit(code=1)
            
            # Try to read all inputs using the registry
            case_inputs = registry.read_case_inputs(case_dir, name="PimpleFoamInputs")
            typer.echo(typer.style("✅ All inputs are valid!", fg=typer.colors.GREEN))
            
            # Show loaded files
            typer.echo("📁 Loaded files:")
            for key, (_, file_spec) in registry.items():
                typer.echo(f"  • {key}: {file_spec.relpath}")

        except ValidationError as e:
            typer.echo(
                typer.style(
                    "❌ Input validation failed:", fg=typer.colors.RED, bold=True
                )
            )
            typer.echo(format_pydantic_errors(e))
            typer.echo(
                typer.style(
                    "Fix the above errors before running the solver.",
                    fg=typer.colors.YELLOW,
                )
            )
            raise typer.Exit(code=1)
        except Exception as e:
            typer.echo(
                typer.style(
                    f"❌ Error validating case: {e}", fg=typer.colors.RED, bold=True
                )
            )
            raise typer.Exit(code=1)

        return

    pimplefoam = PimpleFoam(argv)
    pimplefoam.run()
