#Install PyPDF2
#pip install PyPDF2
#This script takes the first page of 4 PDF files and merges them into one PDF file

import PyPDF2

# specify the PDF file names
pdf1 = "file1.pdf"
pdf2 = "file2.pdf"
pdf3 = "file3.pdf"
pdf4 = "file4.pdf"

# create a pdf file merger object
merger = PyPDF2.PdfMerger()

# open each pdf file
file1 = open(pdf1, "rb")
file2 = open(pdf2, "rb")
file3 = open(pdf3, "rb")
file4 = open(pdf4, "rb")

# append the first page of each pdf file to the merger. (this only merges the first page of each pdf file)
merger.append(fileobj=file1, pages=(0, 1))
merger.merge(position=2, fileobj=file2, pages=(0, 1))
merger.merge(position=3, fileobj=file3, pages=(0, 1))
merger.merge(position=4, fileobj=file4, pages=(0, 1))

# close the input pdf files
file1.close()
file2.close()
file3.close()
file4.close()

# write to an output pdf file
output = open("combined_pages.pdf", "wb")
merger.write(output)
output.close()
