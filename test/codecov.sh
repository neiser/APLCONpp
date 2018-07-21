#!/bin/bash -e

echo "Running build_and_test..."

mkdir -p codecov-build
pushd $(pwd) >/dev/null
cd codecov-build
cmake -DCOVERAGE=On ../.. >/dev/null
nice make build_and_test -j >/dev/null
popd >/dev/null

echo "Running gcovr..."

rm -rf codecov-html
mkdir -p codecov-html
pushd $(pwd) >/dev/null
cd codecov-html
gcovr --object-directory ../codecov-build -r ../.. --html --html-details -o index.html
popd >/dev/null

echo "Done. Open $(pwd)/codecov-html/index.html"
