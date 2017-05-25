# Gaussian-fitting
Matlab codes for Gaussian fitting. ConsGaussFit05082016 is the main function. Others are supporting functions.
Instructions are in the comments of ConsGaussFit05082016.

The format of the input data table: should have one string column of protein ID and a header row. Other columns are protein intensities.

To use the code, set the folder where the codes are put together, to be the working directory.
use command:
ConsGaussFit05082016('SEC_Bio1_nov.csv','sec1-2015nov20160508','plot-uma-2015nov20160508','peakloc-2015nov20160508')

The first argument is input data
Second is output folder name
Third is plot names
Forth is the table of peak locations

The code uses append_pdfs.m and ghostscript.m to combine plots.
The link is:
https://www.mathworks.com/matlabcentral/fileexchange/31215-append-pdfs
