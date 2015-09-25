from openpyxl import load_workbook

class xlsx:
    """
        Class xls

    """

    def __init__(self, xls_file):
        self.path = xls_file
        self.wb = load_workbook(filename=xls_file, read_only=True)
