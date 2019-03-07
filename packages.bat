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


:: install the 'python' ----------------------
ECHO **
ECHO **
ECHO ** install the 'python'
CMD /C " "%SRC_HOME%/nuget.exe"  install python -Version 3.7.1 -OutputDirectory "%SRC_HOME%" "
ECHO **
ECHO **
ECHO ** install the 'pip'
CMD /C " "%SRC_HOME%/python.3.7.1/tools/python" "%SRC_HOME%/get-pip.py"  --no-warn-script-location "
ECHO **
ECHO **
ECHO ** install the 'python packages'
CMD /C " "%SRC_HOME%/python.3.7.1/tools/Scripts/pip.exe" install numpy matplotlib scipy pandas pprint multiprocess times more-itertools concurrent-utils --no-warn-script-location "


GOTO :EndProcess



:: wait to Enter => Good installation
:EndProcess
  SET /P DUMMY=End of installation. Hit ENTER to continue...

