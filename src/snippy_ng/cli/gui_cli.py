import click


def _missing_gradio_message(missing_name: str | None = None) -> str:
    detail = f" Missing Python module: {missing_name}.\n\n" if missing_name else "\n\n"
    return (
        "The Snippy-NG GUI requires the optional Gradio dependency."
        f"{detail}"
        "Install it with:\n"
        "  pip install 'snippy-nextgen[gui]'\n\n"
        "or add Gradio to an existing environment with:\n"
        "  pip install gradio"
    )

@click.command(context_settings={"show_default": True})
@click.option("--host", default="127.0.0.1", help="Host interface for the GUI server.")
@click.option("--port", default=7860, type=click.INT, help="Port for the GUI server.")
@click.option("--share", is_flag=True, default=False, help="Create a public Gradio share link.")
@click.option("--no-browser", is_flag=True, default=False, help="Do not open the GUI in a browser.")
@click.option("--temp-output", is_flag=True, default=False, help="Write GUI run outputs under the Gradio/system temp directory.")
@click.option("--max-cpus", type=click.IntRange(min=1), default=None, help="Maximum CPU value allowed in the GUI.")
@click.option("--no-server-paths", is_flag=True, default=False, help="Hide server-side path inputs and allow uploads only.")
def gui(host: str, port: int, share: bool, no_browser: bool, temp_output: bool, max_cpus: int | None, no_server_paths: bool):
    """
    Launch the Snippy-NG graphical interface.
    """
    try:
        import gradio as gr  # noqa: F401
    except ModuleNotFoundError as exc:
        click.echo(_missing_gradio_message(exc.name), err=True)
        raise SystemExit(1)
    from snippy_ng.gui import create_app
    
    app = create_app(temp_output=temp_output, server_paths=not no_server_paths, max_cpus=max_cpus)

    app.launch(
        server_name=host,
        server_port=port,
        share=share,
        inbrowser=not no_browser,
    )
