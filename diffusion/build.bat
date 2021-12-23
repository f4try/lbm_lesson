@echo off
cd %~dp0

@REM @REM Build for Visual Studio compiler. Run your copy of amd64/vcvars32.bat to setup 64-bit command-line compiler.

@set SOURCES=diffusion.cpp
@set OUT_EXE=diffusion_cpp.exe
@set CFLAGS=-std=c++17 -march=native -O2 -pipe

@REM @set OUT_DIR=Debug
@REM IF NOT EXIST %OUT_DIR% mkdir %OUT_DIR%
@REM @REM cl /nologo /std:c++17 /EHsc /nologo /Zi /MD  %SOURCES% /Fe%OUT_DIR%/%OUT_EXE%.exe /Fo%OUT_DIR%/ /link %LIBS%


@set OUT_DIR=Release
IF NOT EXIST %OUT_DIR% mkdir %OUT_DIR%
clang-cl /nologo /EHsc /Zi /MD /Ox /Oi %SOURCES% /Fe%OUT_DIR%/%OUT_EXE% /Fo%OUT_DIR%/ 
@REM g++ %CFLAGS% %SOURCES% -o %OUT_DIR%/%OUT_EXE% -lm
@REM clang++ %CFLAGS% %SOURCES% -o %OUT_DIR%/%OUT_EXE%
%~dp0/%OUT_DIR%/%OUT_EXE%