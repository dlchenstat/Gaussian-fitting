# Matlab code for Gaussian Fitting in mass-spectrometry data

This is the code for the paper
McBride, Z., Chen, D., Reick, C., Xie, J., & Szymanski, D. B. (2017). Global analysis of membrane-associated protein oligomerization using protein correlation profiling. Molecular & Cellular Proteomics, mcp–000276
http://www.mcponline.org/content/early/2017/09/08/mcp.RA117.000276.full.pdf

This Matlab code is modified from Kristensen (2012).

Matlab codes for Gaussian fitting. ConsGaussFit05082016 is the main function. Others are supporting functions.
Instructions are in the comments of ConsGaussFit05082016.

# Data Format
The format of the input data table: should have one string column of protein ID and a header row. Other columns are protein intensities.

# Instructions
To use the code, set the folder where the codes are put together, to be the working directory.
use command:
ConsGaussFit05082016('SEC_Bio1_nov.csv','sec1-2015nov20160508','plot-2015nov20160508','peakloc-2015nov20160508')

The first argument is input data
Second is output folder name
Third is plot names
Forth is the table of peak locations

# Description of files
consExp.m, consExp3.m, consExp4.m are the functions of 2, 3, 4 Gaussian peaks. Peaks are separated 

normprob.m is used to calculate peak areas.

Please download all code files from this page.
The link is:
https://www.mathworks.com/matlabcentral/fileexchange/31215-append-pdfs

The code contains peakfinder.m, which is another way to find peaks. Also please download all code files from this page.
https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder-x0--sel--thresh--extrema--includeendpoints--interpolate-?

# References

Kristensen, A. R., Gsponer, J., & Foster, L. J. (2012). A high-throughput approach for measuring temporal changes in the interactome. Nature methods, 9(9), 907-909.
