"""
"""
from data import EXPORTED_PING_SPREADSHEET
from ping.apps.upload import PINGUploadSession

if __name__ == '__main__':
    sess = PINGUploadSession()
    sess.login()
    sess.upload_spreadsheet(EXPORTED_PING_SPREADSHEET)
