
import os
import glob
import csv
import xlsxwriter

def csv2xlsx(csvfile):
    print 'Converting %s' % csvfile
    output_filename = os.path.splitext(csvfile)[0] + '.xlsx'
    wb = xlsxwriter.Workbook(output_filename, {'constant_memory': True})
    ws = wb.add_worksheet()
    with open(csvfile, 'rb') as f:
        reader = csv.reader(f)
        for i, row in enumerate(reader):
            ws.write_row(i, 0, row)
    print 'Saving ...'
    wb.close()

def main():
    folders = [x[0] for x in os.walk('.') if x[0] != '.']
    for folder in folders:
        csvs = glob.glob(os.path.join(folder, '*.csv'))
        for csvfile in csvs:
            csv2xlsx(csvfile)

main()
