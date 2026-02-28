def define_env(env):
    @env.macro
    def plotly_embed(src, id, height="775px"):
        """
        This is a macro function added by Claude to allow the embedding of plotly plots which can be expanded to full screen.

        Embed a Plotly HTML file with a fullscreen button.

        Parameters
        ----------
        src : str
            Path to the Plotly HTML file (relative to the page).
        id : str
            A unique HTML id for the iframe (must be unique per page).
        height : str
            CSS height of the embedded iframe, e.g. "750px" or "60vh".
        """
        return (
            f'<div style="display:flex; justify-content:flex-end; margin:.25rem 0;">\n'
            f"<button onclick=\"document.getElementById('{id}').requestFullscreen()\"\n"
            f'  style="cursor:pointer; background:none; border:1px solid #ccc; padding:4px 12px; border-radius:4px;">\n'
            f"  Fullscreen Plot &#x2922;\n"
            f"</button>\n"
            f"</div>\n"
            f"<iframe\n"
            f'  id="{id}"\n'
            f'  src="{src}"\n'
            f'  style="width:100%; height:{height}; border:none;"\n'
            f"  allowfullscreen></iframe>"
        )
