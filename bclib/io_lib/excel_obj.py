import numpy as np
from openpyxl import load_workbook

class xlsx:
    """
        Class xls

    """

    def __init__(self, xls_file):
        self.path = xls_file
        self.wb = load_workbook(filename=xls_file, read_only=True,data_only=True)
    
    def read_spreadsheet_allrow(self,worsheet,y_range):
                
        ws = self.wb[worsheet]
        n = len(ws.columns[1])
        a = np.zeros((n,len(y_range))) 
                         
        for x in range(2,n):
            for y in y_range:
                print(x,y)
                a[x-2][y-y_range[0]] = ws.cell(row = x , column = y).value
                
         
        return a   
    
    def read_spreadsheet_allcols(self,worsheet,x_range):
                
        ws = self.wb[worsheet]
        n = len(ws.rows[1])
        a = np.zeros((len(x_range),n)) 
                         
        for x in x_range:
            for y in range(2,n):
                print(x,y)
                a[x-x_range[0]][y-1] = ws.cell(row = x , column = y).value
                
         
        return a   