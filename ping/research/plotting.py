"""
"""

import os
import os.path as op


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
                fig_name = getattr(plt.figure(fi), 'name', 'figure%d' % fi)
                fig_path = op.join(output_dir, '%s.png' % fig_name)
                if not op.exists(op.dirname(fig_path)):
                    os.makedirs(op.dirname(fig_path))
                plt.savefig(fig_path)
                plt.close()

    elif plotengine in ['bokeh', 'bokeh-silent']:
        import bokeh.plotting
        import tempfile
        output_dir = output_dir or tempfile.mkdtemp()
        output_name = getattr(ax, 'name', ax.title).replace(':', '-')
        output_file = op.join(output_dir, '%s.html' % output_name)

        if not op.exists(output_dir):
            os.makedirs(output_dir)
        if op.exists(output_file):
            os.remove(output_file)
        bokeh.plotting.output_file(output_file, title=ax.title, mode='inline')
        if plotengine == 'bokeh':
            bokeh.plotting.show(ax)
