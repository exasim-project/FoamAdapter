import typer


app = typer.Typer()

@app.callback()
def main():
    """FoamAdapter CLI"""
    pass

@app.command()
def hello(name: str = "World"):
    """Say hello."""
    typer.echo(f"Hello {name}!")
