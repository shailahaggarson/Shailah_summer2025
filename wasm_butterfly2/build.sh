#!/bin/bash

# Build script for Hofstadter Butterfly WebAssembly

echo "Building Hofstadter Butterfly WebAssembly..."

# Check if emscripten is installed
if ! command -v emcc &> /dev/null; then
    echo "Error: Emscripten not found. Please install Emscripten SDK first."
    echo "Visit: https://emscripten.org/docs/getting_started/downloads.html"
    exit 1
fi

# Compile C++ to WebAssembly
emcc hofstadter.cpp \
    -o hofstadter.js \
    -s WASM=1 \
    -s EXPORTED_FUNCTIONS='["_computeSpectrum", "_freeMemory", "_malloc", "_free"]' \
    -s EXPORTED_RUNTIME_METHODS='["ccall", "cwrap", "getValue", "setValue"]' \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s MODULARIZE=1 \
    -s EXPORT_NAME='HofstadterModule' \
    -O3 \
    -std=c++17

if [ $? -eq 0 ]; then
    echo "Build successful! Generated hofstadter.js and hofstadter.wasm"
    echo "You can now open index.html in a web browser (served from a local server)"
else
    echo "Build failed!"
    exit 1
fi