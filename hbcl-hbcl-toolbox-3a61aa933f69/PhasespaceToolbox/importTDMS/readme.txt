MATLAB TDMS IMPORT UTILITY
Author: Michael Corbett, AFRL/RZPE, michael.corbett@wpafb.af.mil, 937-255-8977
Approved for Public Release on 13 JAN 2011 (Case # 88ABW-2011-0105)

Installation:
    Expand the NI DIAdem Connectivity Library with MATLAB TDM/TDMS example ("Reading TDM/TDMS Files with The MathWorks, Inc. MATLAB® Software", http://zone.ni.com/devzone/cda/epd/p/id/5957) into a folder anywhere on the computer. Maintain the directory structure of the "dev" folder as it is (the rest can be discarded if desired) and add that "dev" folder to the MATLAB path with subfolders. Make sure that importTDMS.m file is also somewhere in the MATLAB path. It should automatically be able to locate the required NI .h and .dll files as long as they are on the MATLAB path. 

Usage:
    General usage of importTDMS can be seen by typing 'help importTDMS' at the MATLAB prompt. The .m file is also documented inline to assist in development or troubleshooting. The utility has not been tested with an extremely diverse set of TDMS files, so bugs may still exist. The NI DLL should be capable of processing TDM files as well as TDMS files, but the importTDMS.m script has not been tested with any TDM files (though it was written such that TDM files can be selected in the file dialog box). 
    The plotTorDAC.m file is an example of a script that makes use of the imported TDMS data. It will not be useful in general since it is designed to process a specific data file structure, however it does provide some examples of how to "search" for things in the MATLAB structure that results from the importTDMS function.