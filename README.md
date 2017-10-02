# Matlab code for Gaussian Fitting in mass-spectrometry data
Matlab codes for Gaussian fitting. ConsGaussFit05082016 is the main function. Others are supporting functions.
Instructions are in the comments of ConsGaussFit05082016.

The format of the input data table: should have one string column of protein ID and a header row. Other columns are protein intensities.

To use the code, set the folder where the codes are put together, to be the working directory.
use command:
ConsGaussFit05082016('SEC_Bio1_nov.csv','sec1-2015nov20160508','plot-2015nov20160508','peakloc-2015nov20160508')

The first argument is input data
Second is output folder name
Third is plot names
Forth is the table of peak locations

The code uses append_pdfs.m and ghostscript.m to combine plots.
The link is:
https://www.mathworks.com/matlabcentral/fileexchange/31215-append-pdfs

The code contains peakfinder.m, which is another way to find peaks.
https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder-x0--sel--thresh--extrema--includeendpoints--interpolate-?
