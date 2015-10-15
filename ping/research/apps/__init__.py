"""
"""

from ping.apps import PINGArgParser


class ResearchArgParser(PINGArgParser):
    """
    """
    axis_choices = ['AI:mean', 'AI:std',
                    'LH_PLUS_RH:mean', 'LH_PLUS_RH:std']

    def add_common_parser_args(self, arglist):
        for arg in arglist:
            if arg in ['prefix']:
                self.add_argument('prefix',
                                  help="prefix to include in the analysis")

            elif arg in ['prefixes']:
                self.add_argument('prefixes',
                                  help="comma-separated list of prefixes to"
                                  " include in the analysis")

            elif arg in ['key']:
                self.add_argument('key', choices=self.axis_choices)

            else:
                super(ResearchArgParser, self).add_common_parser_args([arg])

    def parse_args(self, *args, **kwargs):
        outvals = super(ResearchArgParser, self).parse_args(*args, **kwargs)
        if getattr(outvals, 'prefixes', None):
            outvals.prefixes = outvals.prefixes.split(',')
        return outvals
