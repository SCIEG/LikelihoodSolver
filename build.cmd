@ECHO OFF

rmdir /S /Q target
xcopy /E /I additions target
cd src
mingw32-make
cd ..
