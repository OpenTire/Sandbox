@echo off

set MATLAB=C:\Program Files\MATLAB\R2020a

cd .

chcp 1252

if "%1"=="" ("C:\PROGRA~1\MATLAB\R2020a\bin\win64\gmake"  -f mfeval_wrapper_rtw.mk all) else ("C:\PROGRA~1\MATLAB\R2020a\bin\win64\gmake"  -f mfeval_wrapper_rtw.mk %1)
@if errorlevel 1 goto error_exit

exit 0

:error_exit
echo The make command returned an error of %errorlevel%
exit 1
