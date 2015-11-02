"""
"""

import os


def show_plots(plotengine, ax=None, output_dir=None):
    if plotengine == 'mpld3':
        import mpld3
        mpld3.show()

    elif plotengine == 'matplotlib':
        import matplotlib.pyplot as plt
        if not output_dir:  # None or ""
            plt.show()
        else:
            for fi in plt.get_fignums():
                plt.figure(fi)
                plt.savefig('%s.png' % (getattr(plt.figure(fi), 'name', 'figure%d' % fi)))
                plt.close()

    elif plotengine in ['bokeh', 'bokeh-silent']:
        import bokeh.plotting
        import tempfile
        output_dir = output_dir or tempfile.mkdtemp()
        output_file = os.path.join(output_dir, '%s.html' % getattr(ax, 'name', ax.title))

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        if os.path.exists(output_file):
            os.remove(output_file)
        bokeh.plotting.output_file(output_file, title=ax.title, mode='inline')
        if plotengine == 'bokeh':
            bokeh.plotting.show(ax)
