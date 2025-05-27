#!/bin/bash
# Compile C++ to WebAssembly
emcc butterfly.cpp \
  -std=c++11 \
  -O3 \
  -s WASM=1 \
  -s MODULARIZE=1 \
  -s EXPORT_ES6=1 \
  -s USE_ES6_IMPORT_META=0 \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s ENVIRONMENT=web \
  -s EXPORTED_FUNCTIONS="['_calculateSpectrum', '_gcd']" \
  -s EXPORTED_RUNTIME_METHODS='["ccall", "cwrap"]' \
  -o butterfly.js