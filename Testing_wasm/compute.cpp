#include <emscripten.h>  // Emscripten header

// Function to process data (e.g., add two numbers)
extern "C" {
    EMSCRIPTEN_KEEPALIVE  // Expose this function to JS
    int add(int a, int b) {
        return a + b;
    }
}