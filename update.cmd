@echo off
setlocal
set HERE=%~dp0

echo Updating chemdraft from git...
git -C "%HERE%." pull --ff-only
if errorlevel 1 (
    echo.
    echo Update failed. If you have local edits, resolve them first.
    pause
    exit /b 1
)

echo Reinstalling, in case dependencies changed...
python -m pip install --quiet -e "%HERE%."
if errorlevel 1 (
    echo.
    echo Reinstall failed.
    pause
    exit /b 1
)

set SKILLDIR=%USERPROFILE%\.config\opencode\skills\chemdraft
if not exist "%SKILLDIR%" mkdir "%SKILLDIR%"
copy /y "%HERE%skills\chemdraft\SKILL.md" "%SKILLDIR%\SKILL.md" >nul

if not "%~1"=="" (
    if not exist "%~1\.opencode\skill\chemdraft" mkdir "%~1\.opencode\skill\chemdraft"
    copy /y "%HERE%skills\chemdraft\SKILL.md" "%~1\.opencode\skill\chemdraft\SKILL.md" >nul
    echo Skill updated in %~1
)

echo.
echo Updated. Version:
python -c "import chemdraft; print(chemdraft.__version__)"
pause
