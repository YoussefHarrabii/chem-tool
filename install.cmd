@echo off
setlocal
set HERE=%~dp0

echo ============================================
echo  chemdraft installer
echo ============================================
echo.

where python >nul 2>&1
if errorlevel 1 (
    echo Python was not found on PATH.
    echo Install Python 3.10 or newer from python.org and tick
    echo "Add python.exe to PATH" during setup, then run this again.
    echo.
    pause
    exit /b 1
)

python -c "import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)"
if errorlevel 1 (
    echo Python 3.10 or newer is required. You have:
    python --version
    pause
    exit /b 1
)

echo Installing dependencies. RDKit is about 35 MB, this takes a minute...
python -m pip install --quiet --upgrade rdkit pymupdf pillow
if errorlevel 1 (
    echo.
    echo Dependency install failed. Check your internet connection.
    pause
    exit /b 1
)

echo Installing chemdraft...
python -m pip install --quiet -e "%HERE%."
if errorlevel 1 (
    echo.
    echo Package install failed.
    pause
    exit /b 1
)

set SKILLDIR=%USERPROFILE%\.config\opencode\skills\chemdraft
if not exist "%SKILLDIR%" mkdir "%SKILLDIR%"
copy /y "%HERE%skills\chemdraft\SKILL.md" "%SKILLDIR%\SKILL.md" >nul
echo Skill installed for every opencode project.

if not "%~1"=="" (
    if not exist "%~1\.opencode\skill\chemdraft" mkdir "%~1\.opencode\skill\chemdraft"
    copy /y "%HERE%skills\chemdraft\SKILL.md" "%~1\.opencode\skill\chemdraft\SKILL.md" >nul
    echo Skill also installed into %~1
)

echo.
echo Checking that it draws something...
chemdraft smiles "CCO" -o "%TEMP%\chemdraft-check.png"
if errorlevel 1 (
    echo.
    echo The check failed. Send the text above for diagnosis.
    pause
    exit /b 1
)

echo.
echo Done. A test picture is at %TEMP%\chemdraft-check.png
echo In a chat, ask for a molecule and the agent will call chemdraft.
echo.
pause
