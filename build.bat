cd %~dp0

@REM @REM Build for Visual Studio compiler. Run your copy of amd64/vcvars32.bat to setup 64-bit command-line compiler.

@set SOURCES=flow.cpp
@set OUT_EXE=flow_cpp.exe

@REM @set OUT_DIR=Debug
@REM IF NOT EXIST %OUT_DIR% mkdir %OUT_DIR%
@REM @REM cl /nologo /std:c++17 /EHsc /nologo /Zi /MD  %SOURCES% /Fe%OUT_DIR%/%OUT_EXE%.exe /Fo%OUT_DIR%/ /link %LIBS%


@set OUT_DIR=Release
IF NOT EXIST %OUT_DIR% mkdir %OUT_DIR%
cl /nologo /EHsc /Zi /MD /Ox /Oi %SOURCES% /Fe%OUT_DIR%/%OUT_EXE% /Fo%OUT_DIR%/ 

%~dp0/%OUT_DIR%/%OUT_EXE%