#!/bin/sh
rm -rf build
find . \( -name '*.pyc' -o -name '*~' \) -exec rm -rf {} \;
