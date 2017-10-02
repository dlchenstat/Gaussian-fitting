# Matlab code for Gaussian Fitting in mass-spectrometry data

This is the code for the paper
McBride, Z., Chen, D., Reick, C., Xie, J., & Szymanski, D. B. (2017). Global analysis of membrane-associated protein oligomerization using protein correlation profiling. *Molecular & Cellular Proteomics*, mcp–000276
http://www.mcponline.org/content/early/2017/09/08/mcp.RA117.000276.full.pdf

This Matlab code is adapted from Kristensen et al (2012).

Matlab codes for Gaussian fitting. ConsGaussFit05082016 is the main function. Others are supporting functions.
Instructions are in the comments of ConsGaussFit05082016.

# Description of files
consExp.m, consExp3.m, consExp4.m are the functions of 2, 3, 4 Gaussian peaks. Peaks are separated by ≥ 4 fractions

normprob.m is used to calculate peak areas.

# Data Format of Input Dataset
The format of the input data table: should have one or more string columns of protein ID (can include protein annotations) and a header row. Other columns are protein intensities.

# Instructions
To use the code, set the folder where the codes are put together, to be the working directory.
use command:
```Matlab
ConsGaussFit05082016('SEC_Bio1_nov.csv','sec1-2015nov20160508','plot-2015nov20160508','peakloc-2015nov20160508')
```

The first argument is input data
Second is output folder name
Third is plot names
Forth is the table of peak locations

Please download all code files from this page and put them in the working directory.
The link is:
https://www.mathworks.com/matlabcentral/fileexchange/31215-append-pdfs

The code contains peakfinder.m, which is another way to find peaks. Also please download all code files from this page.
https://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder-x0--sel--thresh--extrema--includeendpoints--interpolate-

# Data Format of Results
The output contains a .csv file of fitted peak locations, heights widths. It also contains plots of fitted curves and individual curves if the protein has multiple peaks.

# References
Kristensen, A. R., Gsponer, J., & Foster, L. J. (2012). A high-throughput approach for measuring temporal changes in the interactome. Nature methods, 9(9), 907-909.
