@echo off
chcp 65001 >nul
set PYTHONIOENCODING=utf-8
cd /d "%~dp0"
rem Downloads ONE file (E15N10, release 2020-2025, about 8 MB) into tools\egms_download\check\ for a structural check.
".venv\Scripts\python.exe" egms_dl.py fetch E15N10 2020-2025
echo.
pause
