@echo off
rem Talanai exhibit demo for Windows. Just run: demo
rem Four beats, no network and no engine needed. See tools/demo.py.
rem Uses TalanaiDock's bundled interpreter when TALANAI_PYTHON points at it.
if defined TALANAI_PYTHON (
  "%TALANAI_PYTHON%" "%~dp0tools\demo.py" %*
) else (
  python "%~dp0tools\demo.py" %*
)
