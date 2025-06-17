@ECHO OFF

pushd %~dp0

REM Command file for Sphinx documentation

if "%SPHINXBUILD%" == "" (
	set SPHINXBUILD=sphinx-build
)
set SOURCEDIR=.
set BUILDDIR=_build
set TARGET=%1

if "%TARGET%"=="" (
	set TARGET="help"
)
if "%TARGET%"=="cleanhtml" (
	set TARGET=html
	set SPHINXOPTS=-E -a
)

REM %SPHINXBUILD% >NUL 2>NUL
REM if errorlevel 9009 (
REM 	echo.
REM 	echo.The 'sphinx-build' command was not found. Make sure you have Sphinx
REM 	echo.installed, then set the SPHINXBUILD environment variable to point
REM 	echo.to the full path of the 'sphinx-build' executable. Alternatively you
REM 	echo.may add the Sphinx directory to PATH.
REM 	echo.
REM 	echo.If you don't have Sphinx installed, grab it from
REM 	echo.http://sphinx-doc.org/
REM 	exit /b 1
REM )

set TORUN=%SPHINXBUILD% -M %TARGET% %SOURCEDIR% %BUILDDIR% %SPHINXOPTS% %O%
echo %TORUN%
%TORUN%
popd
