compile:
emcc compute.cpp -s WASM=1 -o compute.js -s EXPORTED_FUNCTIONS="['_add']" -s EXPORTED_RUNTIME_METHODS="['ccall', 'cwrap']"