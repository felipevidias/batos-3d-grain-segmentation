@echo off
rem Phase 5g full run (decided 2026-11-10): frozen v3, cores 160^3 + halo 20 (200^3 padded tiles),
rem one scale per scan (central r_granulo = 10), GT-free specimen-touching tiles (filled first), <= 3 workers, 11 GB, orphan-rescue stitch. Restartable: re-running
rem skips finished tiles. Progress: Get-Content C:\BA_TOS_DATA\runs\fullvolume\full_c160.log -Wait -Tail 20
setlocal
cd /d C:\Users\USER\batos-fullvolume
set PY=C:\BA_TOS_DATA\envs\batos-py312\Scripts\python.exe
set LOG=C:\BA_TOS_DATA\runs\fullvolume\full_c160.log
echo [%date% %time%] PIPELINE START >> "%LOG%"
"%PY%" -I scripts\run_fullvolume.py run --name full_c160 --block specimen --core 160 --scale-rule per_scan --workers 3 --mem-budget-gb 11 >> "%LOG%" 2>&1 || goto :fail
"%PY%" -I scripts\run_fullvolume.py stitch --name full_c160 --orphan-rescue --variant _rescue >> "%LOG%" 2>&1 || goto :fail
"%PY%" -I scripts\run_fullvolume.py evaluate --name full_c160 --variant _rescue >> "%LOG%" 2>&1 || goto :fail
"%PY%" -I scripts\run_fullvolume.py evaluate --name full_c160 --variant _rescue --gt-set listed >> "%LOG%" 2>&1 || goto :fail
"%PY%" -I scripts\report_fullvolume.py --name full_c160 --variant _rescue >> "%LOG%" 2>&1 || goto :fail
"%PY%" -I scripts\report_fullvolume.py --name full_c160 --variant _rescue --gt-set listed >> "%LOG%" 2>&1 || goto :fail
echo [%date% %time%] PIPELINE DONE >> "%LOG%"
exit /b 0
:fail
echo [%date% %time%] PIPELINE FAILED errorlevel %errorlevel% >> "%LOG%"
exit /b 1
