const sketch = (p) => {
  let addFunction = null;

  p.setup = () => {
    p.createCanvas(400, 400);
    p.textSize(32);
    
    // Get function once ready
    Module.onRuntimeInitialized = () => {
      addFunction = Module.cwrap('add', 'number', ['number', 'number']);
      console.log("3 + 5 =", addFunction(3, 5)); // Test call
    };
  };

  p.draw = () => {
    p.background(220);
    if (addFunction) {
      const result = addFunction(10, 20);
      p.text(`10 + 20 = ${result}`, 50, 100);
    } else {
      p.text("Loading WASM...", 50, 100);
    }
  };
};

new p5(sketch);