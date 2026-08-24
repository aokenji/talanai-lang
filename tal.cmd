@echo off
rem Talanai launcher for Windows. "tal check study.tal"
rem Uses TalanaiDock's bundled interpreter when TALANAI_PYTHON points at it.
if defined TALANAI_PYTHON (
  "%TALANAI_PYTHON%" "%~dp0tal.py" %*
) else (
  python "%~dp0tal.py" %*
)
