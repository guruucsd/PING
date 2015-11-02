"""
Export all measures to a spreadsheet, then
upload to the PING server.
"""
import os
from functools import partial

from .export import EXPORTED_PING_SPREADSHEET
from ..ping.apps import PINGSession
from ..research.data import get_all_data


def do_upload(csv_file=EXPORTED_PING_SPREADSHEET,
              atlas='desikan', data_dir='data',
              username=None, passwd=None, force=False):

    if force or not os.path.exists(csv_file):
        # Limit to imaging data only
        data = get_all_data(atlas, data_dir=data_dir,
                            username=username, passwd=passwd)
        # data.filter(filter_fns=[partial(lambda k, v, p: k.startswith(p), p=p)
        #                         for p in data.IMAGING_PREFIX])
        data.export(out_file=csv_file)

    # Upload the new spreadsheet.
    sess = PINGSession()
    sess.login()
    sess.upload_user_spreadsheet(csv_file)


if __name__ == '__main__':
    parser = ResearchArgParser(description="Upload data spreadsheet to the PING portal.",
                               common_args=['atlas', 'data-dir',
                                            'username', 'passwd', 'force'])
    parser.add_argument('--csv_file', nargs='?', default=EXPORTED_PING_SPREADSHEET)

    args = parser.parse_args()
    do_upload(**vars(args))
