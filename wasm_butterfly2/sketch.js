// Hofstadter Butterfly Visualization using p5.js and WebAssembly

let wasmModule;
let spectrumData = [];
let isLoaded = false;
let isComputing = false;

// Visualization parameters
let maxQ = 50;  // Maximum denominator
let lambda = 2.0;  // Coupling strength
let canvasWidth = 800;
let canvasHeight = 600;

// Display parameters
let xScale, yScale;
let xOffset = 50;
let yOffset = 50;
let plotWidth, plotHeight;

async function preload() {
    console.log('Loading WebAssembly module...');
    try {
        // Load the WebAssembly module
        wasmModule = await HofstadterModule();
        console.log('WebAssembly module loaded successfully');
        isLoaded = true;
    } catch (error) {
        console.error('Failed to load WebAssembly module:', error);
    }
}

function setup() {
    createCanvas(canvasWidth, canvasHeight);
    
    plotWidth = canvasWidth - 2 * xOffset;
    plotHeight = canvasHeight - 2 * yOffset;
    
    // Set up coordinate scaling
    xScale = plotWidth;  // theta goes from 0 to 1
    yScale = plotHeight / 8.0;  // eigenvalues typically range from -4 to 4
    
    // Create UI controls
    createControls();
    
    if (isLoaded) {
        computeSpectrum();
    }
}

function createControls() {
    // Create control panel
    let controlDiv = createDiv('');
    controlDiv.position(20, canvasHeight + 20);
    controlDiv.style('font-family', 'Arial');
    
    // Max Q slider
    let qLabel = createP('Max Q (complexity): ' + maxQ);
    qLabel.parent(controlDiv);
    let qSlider = createSlider(10, 100, maxQ, 5);
    qSlider.parent(controlDiv);
    qSlider.input(() => {
        maxQ = qSlider.value();
        qLabel.html('Max Q (complexity): ' + maxQ);
        if (isLoaded && !isComputing) {
            computeSpectrum();
        }
    });
    
    // Lambda slider
    let lambdaLabel = createP('Lambda (coupling): ' + lambda.toFixed(2));
    lambdaLabel.parent(controlDiv);
    let lambdaSlider = createSlider(0.1, 4.0, lambda, 0.1);
    lambdaSlider.parent(controlDiv);
    lambdaSlider.input(() => {
        lambda = lambdaSlider.value();
        lambdaLabel.html('Lambda (coupling): ' + lambda.toFixed(2));
        if (isLoaded && !isComputing) {
            computeSpectrum();
        }
    });
    
    // Compute button
    let computeBtn = createButton('Recompute');
    computeBtn.parent(controlDiv);
    computeBtn.mousePressed(() => {
        if (isLoaded && !isComputing) {
            computeSpectrum();
        }
    });
}

async function computeSpectrum() {
    if (!isLoaded || isComputing) return;
    
    isComputing = true;
    console.log(`Computing spectrum with maxQ=${maxQ}, lambda=${lambda}...`);
    
    try {
        // Call the WebAssembly function
        const resultSizePtr = wasmModule._malloc(4);  // int size
        const dataPtr = wasmModule.ccall('computeSpectrum', 'number', 
            ['number', 'number', 'number'], 
            [maxQ, lambda, resultSizePtr]);
        
        const resultSize = wasmModule.getValue(resultSizePtr, 'i32');
        wasmModule._free(resultSizePtr);
        
        // Read the data
        spectrumData = [];
        let idx = 0;
        
        while (idx < resultSize) {
            const theta = wasmModule.getValue(dataPtr + idx * 8, 'double');
            const count = wasmModule.getValue(dataPtr + (idx + 1) * 8, 'double');
            
            const eigenvalues = [];
            for (let i = 0; i < count; i++) {
                const eig = wasmModule.getValue(dataPtr + (idx + 2 + i) * 8, 'double');
                eigenvalues.push(eig);
            }
            
            spectrumData.push({
                theta: theta,
                eigenvalues: eigenvalues
            });
            
            idx += 2 + count;
        }
        
        // Free the allocated memory
        wasmModule.ccall('freeMemory', null, ['number'], [dataPtr]);
        
        console.log(`Computed ${spectrumData.length} spectrum lines`);
        
    } catch (error) {
        console.error('Error computing spectrum:', error);
    }
    
    isComputing = false;
}

function draw() {
    background(0);
    
    if (!isLoaded) {
        // Show loading message
        fill(255);
        textAlign(CENTER, CENTER);
        textSize(20);
        text('Loading WebAssembly module...', width/2, height/2);
        return;
    }
    
    if (isComputing) {
        fill(255);
        textAlign(CENTER, CENTER);
        textSize(16);
        text('Computing spectrum...', width/2, height/2);
        return;
    }
    
    // Draw axes
    drawAxes();
    
    // Draw spectrum data
    drawSpectrum();
    
    // Draw title and info
    drawTitle();
}

function drawAxes() {
    stroke(100);
    strokeWeight(1);
    
    // X-axis
    line(xOffset, height - yOffset, width - xOffset, height - yOffset);
    
    // Y-axis
    line(xOffset, yOffset, xOffset, height - yOffset);
    
    // X-axis labels
    fill(150);
    textAlign(CENTER, TOP);
    textSize(10);
    for (let i = 0; i <= 10; i++) {
        let x = xOffset + (i * plotWidth) / 10;
        let label = (i / 10).toFixed(1);
        text(label, x, height - yOffset + 5);
        
        // Tick marks
        stroke(100);
        line(x, height - yOffset, x, height - yOffset + 3);
    }
    
    // Y-axis labels
    textAlign(RIGHT, CENTER);
    for (let i = -4; i <= 4; i++) {
        let y = height/2 - i * yScale/2;
        if (y >= yOffset && y <= height - yOffset) {
            text(i.toString(), xOffset - 5, y);
            
            // Tick marks
            stroke(100);
            line(xOffset - 3, y, xOffset, y);
        }
    }
}

function drawSpectrum() {
    strokeWeight(1);
    
    for (let data of spectrumData) {
        let x = xOffset + data.theta * plotWidth;
        
        // Use different colors for different density of eigenvalues
        let numEigs = data.eigenvalues.length;
        let hue = map(numEigs, 1, 20, 240, 0);  // Blue to red
        colorMode(HSB);
        stroke(hue, 80, 90, 150);
        
        // Draw vertical line segments for each eigenvalue interval
        for (let i = 0; i < data.eigenvalues.length - 1; i++) {
            let y1 = height/2 - data.eigenvalues[i] * yScale/2;
            let y2 = height/2 - data.eigenvalues[i+1] * yScale/2;
            
            // Clamp to visible area
            y1 = constrain(y1, yOffset, height - yOffset);
            y2 = constrain(y2, yOffset, height - yOffset);
            
            line(x, y1, x, y2);
        }
        
        // Also draw individual eigenvalue points
        stroke(hue, 100, 100, 200);
        for (let eig of data.eigenvalues) {
            let y = height/2 - eig * yScale/2;
            if (y >= yOffset && y <= height - yOffset) {
                point(x, y);
            }
        }
    }
    
    colorMode(RGB);
}

function drawTitle() {
    fill(255);
    textAlign(LEFT, TOP);
    textSize(16);
    text('Hofstadter Butterfly', 20, 20);
    
    textSize(12);
    text(`Max Q: ${maxQ}, λ: ${lambda.toFixed(2)}`, 20, 45);
    text(`Data points: ${spectrumData.length}`, 20, 65);
    
    // Legend
    textAlign(RIGHT, TOP);
    text('θ (flux ratio) →', width - 20, height - 30);
    
    push();
    translate(30, height/2);
    rotate(-PI/2);
    textAlign(CENTER, BOTTOM);
    text('Energy →', 0, 0);
    pop();
}

// Handle window resize
function windowResized() {
    // Keep the canvas size fixed for consistent visualization
}