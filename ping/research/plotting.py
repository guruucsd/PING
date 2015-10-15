"""
"""

def show_plots(plotengine, ax=None):
    if plotengine == 'mpld3':
        import mpld3
        mpld3.show()
    elif plotengine == 'matplotlib':
        import matplotlib.pyplot as plt
        plt.show()
    elif plotengine == 'bokeh':
        import bokeh.plotting
        import tempfile
        bokeh.plotting.output_file(tempfile.mkstemp()[1] + ".html",
                                   title=ax.title, mode='absolute-dev')
        bokeh.plotting.show(ax)
