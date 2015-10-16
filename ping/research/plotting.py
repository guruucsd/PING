"""
"""

import os


def show_plots(plotengine, ax=None, output_dir=None):
    if plotengine == 'mpld3':
        import mpld3
        mpld3.show()

    elif plotengine == 'matplotlib':
        import matplotlib.pyplot as plt
        plt.show()

    elif plotengine == 'bokeh':
        import bokeh.plotting
        import tempfile
        output_dir = output_dir or tempfile.mkdtemp()
        output_file = os.path.join(output_dir, '%s.html' % getattr(ax, 'name', ax.title))

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        if os.path.exists(output_file):
            os.remove(output_file)
        bokeh.plotting.output_file(output_file, title=ax.title, mode='absolute-dev')
        bokeh.plotting.show(ax)
