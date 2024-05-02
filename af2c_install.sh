#!/bin/bash

# Exit when any command fails
set -e

rm -rf af2complex
git clone https://github.com/FreshAirTonight/af2complex.git
cd af2complex
git checkout ca50922a668d278cd92ad8adb843cd5206f8e8fc

patch=../af2complex.patch
git apply "$patch"
