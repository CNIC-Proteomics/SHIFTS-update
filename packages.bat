@ECHO OFF

:: sets the env. variables from input parameters ----------------------
ECHO **
ECHO **
ECHO ** sets the env. variables from input parameters:

ECHO **
SET  SRC_HOME="%~dp0"
SET  SRC_HOME=%SRC_HOME:"=%
SET  SRC_HOME=%SRC_HOME:~0,-1%

ECHO **
SET  LIB_VERSION=0.1
SET  LIB_PATH="%SRC_HOME%"

SET  LIB_HOME=%SRC_HOME%
ECHO %LIB_HOME%


REM :: install the 'python' ----------------------
REM ECHO **
REM ECHO **
REM ECHO ** install the 'python'
CMD /C " "%SRC_HOME%/nuget.exe"  install python -Version 3.7.1 -OutputDirectory "%SRC_HOME%" "
CMD /C " "%SRC_HOME%/python.3.7.1/tools/python" "%SRC_HOME%/get-pip.py"  --no-warn-script-location "
CMD /C " "%SRC_HOME%/python.3.7.1/tools/Scripts/pip.exe" install numpy matplotlib scipy pandas pprint multiprocess times more-itertools concurrent-utils --no-warn-script-location "

SET /P DUMMY=End of installation. Hit ENTER to continue...

REM :: install the PIP package ----------------------
REM ECHO **
REM ECHO **
REM ECHO ** install the 'pip' package
REM CMD /C " "%PYTHON3x_HOME%/tools/python" "%SRC_HOME%/install/get-pip.py"  --no-warn-script-location "


REM :: install required packages ----------------------
REM ECHO **
REM ECHO **
REM ECHO ** install required packages
REM CMD /C " "%PYTHON3x_HOME%/tools/Scripts/pip3.exe" install numpy matplotlib scipy snakemake pandas pprint --no-warn-script-location "


REM :: install R packages ----------------------
REM ECHO **
REM ECHO **
REM ECHO ** install R packages
REM CMD /C " "%R_HOME%/bin/R" --vanilla --args "%R_LIB%" test2=no < "%SRC_HOME%/install/src/install_Rlibs.R" "


REM :: download and install npm ----------------------
REM ECHO **
REM ECHO **
REM ECHO ** download and install npm
REM SET  NODE_URL=https://nodejs.org/dist/v10.14.2/node-v10.14.2-win-x64.zip
REM CMD /C " "%PYTHON3x_HOME%/tools/python" "%SRC_HOME%/install/src/install_url_pkg.py" "%NODE_URL%" "%NODE_HOME%" "%LIB_HOME%/tmp" move "


REM :: install electron package ----------------------
REM ECHO **
REM ECHO **
REM ECHO ** install electron package
REM CMD /C " "%NODE_HOME%/npm" config set scripts-prepend-node-path true"
REM CMD /C " "%NODE_HOME%/npm" install electron --save-dev --save-exact --global "
REM CMD /C " "%NODE_HOME%/npm" install ps-tree --global "




REM GOTO :EndProcess



REM :: wait to Enter => Good installation
REM :EndProcess
    REM SET /P DUMMY=End of installation. Hit ENTER to continue...
REM © 2019 GitHub, Inc.
REM Terms
REM Privacy
REM Security
REM Status
REM Help
REM Contact GitHub
REM Pricing
REM API
REM Training
REM Blog
REM About
