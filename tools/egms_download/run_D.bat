@echo off
chcp 65001 >nul
set PYTHONIOENCODING=utf-8
rem EGMS allows at most 2 parallel downloads per user (2026-11-20 note); keep workers at 2.
set EGMS_WORKERS=2
rem Confirm03 part D downloads (v4). Reads the CLMS service key from env var EGMS_SERVICE_KEY_JSON (never printed or written).
rem Runs up to 12 passes, 10 minutes apart, and stops early once DONE_D.txt exists.
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo [setup] private venv missing; please run run.bat once first, or tell Claude.
  goto end
)
set /a PASS=0
:loop
set /a PASS+=1
echo.
echo ===== pass %PASS% / 12 =====
".venv\Scripts\python.exe" egms_dl_d.py all
if exist "DONE_D.txt" goto end
if %PASS% GEQ 12 goto end
echo Waiting 10 minutes before the next pass (close this window to stop)...
timeout /t 600 /nobreak >nul
goto loop
:end
echo.
echo Finished. DONE_D.txt present = all part-D files downloaded. Otherwise double-click run_D.bat again later.
pause
