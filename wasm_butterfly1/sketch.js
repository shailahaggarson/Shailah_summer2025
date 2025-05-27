let qMax = 50;
let lambda = 2.0;
let isDrawing = false;
let points = [];
let moduleLoaded = false;

function setup() {
  let canvas = createCanvas(800, 600);
  canvas.parent('canvas-container');
  background(255);
  drawAxes();
  
  // Setup event listeners
  document.getElementById('drawBtn').addEventListener('click', () => {
    if (!moduleLoaded) {
      document.getElementById('status').textContent = "WASM module not loaded yet";
      return;
    }
    drawButterfly();
  });
  
  document.getElementById('clearBtn').addEventListener('click', () => {
    background(255);
    drawAxes();
    document.getElementById('status').textContent = "Cleared";
  });
  
  const qMaxSlider = document.getElementById('qMax');
  const qMaxValue = document.getElementById('qMaxValue');
  qMaxSlider.addEventListener('input', () => {
    qMax = parseInt(qMaxSlider.value);
    qMaxValue.textContent = qMax;
  });
  
  const lambdaSlider = document.getElementById('lambda');
  const lambdaValue = document.getElementById('lambdaValue');
  lambdaSlider.addEventListener('input', () => {
    lambda = parseFloat(lambdaSlider.value) / 10;
    lambdaValue.textContent = lambda.toFixed(1);
  });
}

async function drawButterfly() {
  if (isDrawing || !ButterflyCalc) return;
  isDrawing = true;
  const status = document.getElementById('status');
  status.textContent = "Calculating...";
  
  background(255);
  drawAxes();
  
  // Draw the zero line first
  drawSpectrum(0, 1, lambda, color(0, 0, 255));
  
  // Draw spectra for fractions p/q where q <= qMax
  let drawn = 0;
  const total = Array.from({length: qMax}, (_, q) => {
    return Array.from({length: q+1}, (_, p) => p)
      .filter(p => ButterflyCalc.gcd(p, q+1) === 1).length;
  }).reduce((a, b) => a + b, 0);
  
  for (let q = 1; q <= qMax; q++) {
    for (let p = 1; p <= q; p++) {
      if (ButterflyCalc.gcd(p, q) === 1) {
        drawSpectrum(p, q, lambda, color(0, 100));
        drawn++;
        if (drawn % 10 === 0) {
          status.textContent = `Calculating... ${Math.round(drawn/total*100)}%`;
          await new Promise(resolve => setTimeout(resolve, 0));
        }
      }
    }
  }
  
  status.textContent = "Done!";
  isDrawing = false;
}

function drawSpectrum(p, q, lambda, col) {
  // Calculate eigenvalues using WASM
  const spectrum = ButterflyCalc.calculateSpectrum(p, q, lambda);
  const Xper = spectrum[0];
  const Xanti = spectrum[1];
  
  // Map to screen coordinates
  const theta = p / q;
  const x = map(theta, 0, 1, 50, width - 50);
  
  // Draw the spectrum lines
  stroke(col);
  strokeWeight(0.5);
  
  for (let i = 0; i < Xper.length; i++) {
    const y1 = map(Xper[i], -4, 4, height - 50, 50);
    const y2 = map(Xanti[i], -4, 4, height - 50, 50);
    line(x, y1, x, y2);
  }
}

function drawAxes() {
  stroke(0);
  strokeWeight(1);
  // X-axis (magnetic flux)
  line(50, height - 50, width - 50, height - 50);
  // Y-axis (energy)
  line(50, 50, 50, height - 50);
  
  // X-axis labels
  textSize(12);
  textAlign(CENTER, TOP);
  noStroke();
  fill(0);
  for (let i = 0; i <= 10; i++) {
    let x = map(i/10, 0, 1, 50, width - 50);
    text(nf(i/10, 1, 1), x, height - 40);
    stroke(0);
    line(x, height - 50, x, height - 45);
    noStroke();
  }
  
  // Y-axis labels
  textAlign(RIGHT, CENTER);
  for (let i = -4; i <= 4; i++) {
    let y = map(i, -4, 4, height - 50, 50);
    text(i, 45, y);
    stroke(0);
    line(50, y, 55, y);
    noStroke();
  }
  
  // Axis titles
  textAlign(CENTER, CENTER);
  text("Magnetic Flux (θ = p/q)", width/2, height - 20);
  push();
  translate(20, height/2);
  rotate(-HALF_PI);
  text("Energy Spectrum", 0, 0);
  pop();
  
  // Title
  textSize(16);
  text("Hofstadter Butterfly (WASM)", width/2, 30);
}

function draw() {
  // Animation could go here if needed
}