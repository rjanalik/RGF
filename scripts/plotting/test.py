#!/usr/bin/env python3
import subprocess

error = subprocess.run(['./test.sh World'], shell=True, capture_output=True)
print(error.stdout)
print(error.stderr)
