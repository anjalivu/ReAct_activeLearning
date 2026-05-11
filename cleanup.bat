@echo off
title Cleanup - Keep Only Specific Files
color 0A

echo ========================================
echo  SELECTIVE FILE CLEANUP TOOL
echo ========================================
echo.
echo This script will DELETE all files EXCEPT:
echo   - .odb files
echo   - .dat files
echo   - .xls, .xlsx, .xlsm (Excel files)
echo   - .csv files
echo   - .py (Python source files)
echo   - .ipynb (Jupyter notebooks)
echo   - All FOLDERS (preserved)
echo.
echo Current folder: %CD%
echo.
echo ========================================

REM Show what will be kept
echo.
echo Files that will be KEPT:
echo.
for %%F in (*.odb *.dat *.xls *.xlsx *.xlsm *.csv *.py *.ipynb) do (
    if exist "%%F" echo   [KEEP] %%F
)

REM Count folders
echo.
echo Folders (will be preserved):
for /d %%D in (*) do (
    echo   [FOLDER] %%D
)

echo.
echo ========================================
echo.
echo WARNING: All other files will be PERMANENTLY DELETED!
echo.
set /p CONFIRM=Type YES to continue (or anything else to cancel): 

if /I not "%CONFIRM%"=="YES" (
    echo.
    echo Operation CANCELLED.
    echo.
    pause
    exit /b
)

echo.
echo ========================================
echo  DELETING FILES...
echo ========================================
echo.

setlocal enabledelayedexpansion
set COUNT=0

for %%F in (*) do (
    set "DELETE=1"
    
    REM Don't delete the batch file itself
    if /I "%%F"=="%~nx0" set "DELETE=0"
    
    REM Don't delete files with these extensions
    if /I "%%~xF"==".odb" set "DELETE=0"
    if /I "%%~xF"==".dat" set "DELETE=0"
    if /I "%%~xF"==".xls" set "DELETE=0"
    if /I "%%~xF"==".xlsx" set "DELETE=0"
    if /I "%%~xF"==".xlsm" set "DELETE=0"
    if /I "%%~xF"==".csv" set "DELETE=0"
    if /I "%%~xF"==".py" set "DELETE=0"
    if /I "%%~xF"==".ipynb" set "DELETE=0"
    
    REM Delete if not protected
    if "!DELETE!"=="1" (
        del "%%F" 2>nul
        if errorlevel 1 (
            echo [FAILED] %%F
        ) else (
            echo [DELETED] %%F
            set /a COUNT+=1
        )
    )
)

echo.
echo ========================================
echo  CLEANUP COMPLETE!
echo ========================================
echo.
echo Total files deleted: %COUNT%
echo.
pause