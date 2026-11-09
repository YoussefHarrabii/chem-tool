@echo off
setlocal enabledelayedexpansion
title chemdraft setup

set REPO=YoussefHarrabii/chem-tool
set DEST=%LOCALAPPDATA%\chemdraft
set ZIP=%TEMP%\chemdraft-main.zip

echo ============================================
echo   chemdraft - chemistry diagrams
echo ============================================
echo.

rem ---- Python -----------------------------------------------------------
where python >nul 2>&1
if errorlevel 1 (
    echo Python is not installed. Trying to install it automatically...
    winget install -e --id Python.Python.3.12 --scope user ^
        --accept-package-agreements --accept-source-agreements >nul 2>&1
    set PATH=%LOCALAPPDATA%\Programs\Python\Python312;%LOCALAPPDATA%\Programs\Python\Python312\Scripts;%PATH%
)
where python >nul 2>&1
if errorlevel 1 (
    echo.
    echo Python is needed once. Install it from https://www.python.org/downloads/
    echo During setup tick "Add python.exe to PATH", then double-click this
    echo file again.
    echo.
    pause
    exit /b 1
)

rem ---- download ---------------------------------------------------------
echo Downloading chemdraft...
powershell -NoProfile -Command "try { Invoke-WebRequest -UseBasicParsing 'https://github.com/%REPO%/archive/refs/heads/main.zip' -OutFile '%ZIP%' } catch { exit 1 }"
if errorlevel 1 (
    echo.
    echo Could not download. Check the internet connection.
    echo If the repository was made private, this one-file installer will not
    echo work; use the git instructions in SETUP.md instead.
    echo.
    pause
    exit /b 1
)

echo Unpacking...
if exist "%DEST%" rmdir /s /q "%DEST%" >nul 2>&1
powershell -NoProfile -Command "Expand-Archive -Force '%ZIP%' '%DEST%'" >nul 2>&1

set SRC=
for /d %%D in ("%DEST%\*") do set SRC=%%D
if not exist "%SRC%\pyproject.toml" (
    echo.
    echo Unpacking failed.
    pause
    exit /b 1
)

rem ---- install ----------------------------------------------------------
echo Installing. RDKit is about 35 MB, this takes a minute...
python -m pip install --quiet --upgrade rdkit pymupdf pillow
if errorlevel 1 (
    echo.
    echo Could not install the Python packages. Check the internet connection.
    pause
    exit /b 1
)

python -m pip install --quiet -e "%SRC%"
if errorlevel 1 (
    echo.
    echo Could not install chemdraft.
    pause
    exit /b 1
)

set SKILLDIR=%USERPROFILE%\.config\opencode\skills\chemdraft
if not exist "%SKILLDIR%" mkdir "%SKILLDIR%"
copy /y "%SRC%\skills\chemdraft\SKILL.md" "%SKILLDIR%\SKILL.md" >nul

rem ---- check ------------------------------------------------------------
echo.
echo Checking that it works...
chemdraft smiles "CCO" -o "%TEMP%\chemdraft-check.png" >nul 2>&1
if errorlevel 1 (
    echo.
    echo The check failed. Copy the text above and send it back.
    pause
    exit /b 1
)

echo.
echo Done. A test picture is at:
echo   %TEMP%\chemdraft-check.png
echo.
echo In opencode, just ask for a molecule or a reaction.
echo Double-click this same file whenever you want the latest version.
echo.
pause
