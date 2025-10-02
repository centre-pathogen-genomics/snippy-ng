import click 

def error(msg):
    click.echo(f"Error: {msg}", err=True)
    raise click.Abort()