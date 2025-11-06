# Developer Manual - logs.py

These functions track user actions for verification purposes rather than developer debugging. Logs are written both to the terminal (unless *logs.VERBOSE=False*, in which case only warnings are shown) and to a log file located in *~usr/.CALINS/*. This folder is checked at each CALINS import, and the last file is deleted if the folder contains more than 10 files (since a new file is created each day CALINS is imported).

The *log_exec* decorator allows tracking the execution of CALINS's main functions.