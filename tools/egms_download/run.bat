@echo off
chcp 65001 >nul
set PYTHONIOENCODING=utf-8
rem Reads the CLMS service key from env var EGMS_SERVICE_KEY_JSON (JSON text or a path to the JSON file).
rem To use another variable name, remove "rem" below and put the NAME only:
rem set EGMS_VAR=NAME
cd /d "%~dp0"
set BASEPY=C:\Users\Jupyter\anaconda3\envs\llama4_env\python.exe
if not exist "%BASEPY%" set BASEPY=python
if not exist ".venv\Scripts\python.exe" (
  echo [setup] creating private venv in tools\egms_download\.venv
  "%BASEPY%" -m venv .venv
  if errorlevel 1 goto fail
)
".venv\Scripts\python.exe" -c "import jwt, cryptography" 2>nul
if errorlevel 1 (
  echo [setup] installing pyjwt and cryptography into the private venv
  ".venv\Scripts\python.exe" -m pip install --disable-pip-version-check -q pyjwt cryptography
  if errorlevel 1 goto fail
)
".venv\Scripts\python.exe" egms_dl.py all
goto end
:fail
echo [setup] failed. Please send the messages above to Claude.
:end
echo.
echo Finished. If it stopped midway, double-click run.bat again to resume.
pause
