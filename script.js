let chartInstance = null;

function readXYZFile(file, callback) {
    const reader = new FileReader();
    reader.onload = function (e) {
        try {
            const lines = e.target.result.split('\n');
            const numParticles = parseInt(lines[0].trim());
            const boxLine = lines[1].trim().split(/\s+/);
            
            if (boxLine[0] !== "Timestep") {
                throw new Error("Expected 'Timestep' in second line of .xyz file");
            }

            const boxIndex = boxLine.indexOf("box");
            if (boxIndex === -1) {
                throw new Error("'box' keyword not found in second line of .xyz file");
            }

            const boxData = boxLine.slice(boxIndex + 1, boxIndex + 7).map(parseFloat);
            if (boxData.length < 6) {
                throw new Error("Insufficient data for box dimensions after 'box' keyword");
            }

            const boxSize = [
                boxData[1] - boxData[0], // xhi - xlo
                boxData[3] - boxData[2], // yhi - ylo
                boxData[5] - boxData[4]  // zhi - zlo
            ];

            const positions = [];
            for (let i = 2; i < 2 + numParticles; i++) {
                const parts = lines[i].trim().split(/\s+/);
                if (parts.length >= 4) {
                    positions.push(parts.slice(1, 4).map(parseFloat));
                }
            }

            callback(positions, boxSize);
        } catch (error) {
            document.getElementById('error').textContent = error.message;
        }
    };
    reader.readAsText(file);
}

function dot(a, b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function computeADF(positions, boxSize, rMax, numBins) {
    const numParticles = positions.length;
    const thetaBin = Array.from({ length: numBins + 1 }, (_, i) => (i * 180) / numBins);
    const dtheta = thetaBin[1] - thetaBin[0];
    const adf = new Array(numBins).fill(0);

    for (let i = 0; i < numParticles; i++) {
        const delta = positions.map(pos => [
            pos[0] - positions[i][0],
            pos[1] - positions[i][1],
            pos[2] - positions[i][2]
        ]);
        
        // Apply periodic boundary conditions
        delta.forEach((d, idx) => {
            for (let j = 0; j < 3; j++) {
                d[j] = d[j] - boxSize[j] * Math.round(d[j] / boxSize[j]);
            }
        });

        const distances = delta.map(d => Math.sqrt(dot(d, d)));
        const neighbors = distances
            .map((d, idx) => (d > 0 && d < rMax ? idx : -1))
            .filter(idx => idx !== -1);

        for (let jIdx = 0; jIdx < neighbors.length; jIdx++) {
            for (let kIdx = jIdx + 1; kIdx < neighbors.length; kIdx++) {
                const j = neighbors[jIdx];
                const k = neighbors[kIdx];
                const vecIJ = delta[j];
                const vecIK = delta[k];
                const normIJ = Math.sqrt(dot(vecIJ, vecIJ));
                const normIK = Math.sqrt(dot(vecIK, vecIK));
                if (normIJ === 0 || normIK === 0) continue;

                let cosTheta = dot(vecIJ, vecIK) / (normIJ * normIK);
                cosTheta = Math.min(1.0, Math.max(-1.0, cosTheta));
                const theta = (Math.acos(cosTheta) * 180) / Math.PI;
                const binIndex = Math.floor(theta / dtheta);
                if (binIndex >= 0 && binIndex < numBins) {
                    adf[binIndex] += 2;
                }
            }
        }
    }

    const totalCounts = adf.reduce((sum, val) => sum + val, 0);
    if (totalCounts === 0) {
        console.warn("No angles were counted within r_max.");
        return { theta: thetaBin.slice(0, -1), adf };
    }

    const dthetaRad = (dtheta * Math.PI) / 180;
    const normalizedADF = adf.map(count => count / (totalCounts * dthetaRad));
    return { theta: thetaBin.slice(0, -1), adf: normalizedADF };
}

function plotADF(theta, adf) {
    const ctx = document.getElementById('adfChart').getContext('2d');
    
    if (chartInstance) {
        chartInstance.destroy();
    }

    chartInstance = new Chart(ctx, {
        type: 'line',
        data: {
            labels: theta.map(val => val.toFixed(2)),
            datasets: [{
                label: 'ADF',
                data: adf,
                borderColor: '#007bff',
                backgroundColor: 'rgba(0, 123, 255, 0.1)',
                fill: false,
                pointRadius: 0,
                pointHoverRadius: 3
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                legend: { display: true },
                title: {
                    display: true,
                    text: 'Angular Distribution Function'
                }
            },
            scales: {
                x: {
                    title: { display: true, text: 'Angle (degrees)' },
                    grid: { display: true }
                },
                y: {
                    title: { display: true, text: 'Density (1/radian)' },
                    grid: { display: true }
                }
            }
        }
    });
}

// Event listeners
document.getElementById('fileInput').addEventListener('change', function (e) {
    const file = e.target.files[0];
    if (!file) return;
    document.getElementById('plotButton').disabled = false;
});

document.getElementById('plotButton').addEventListener('click', function () {
    const file = document.getElementById('fileInput').files[0];
    const rMax = parseFloat(document.getElementById('rMaxInput').value);
    const numBins = parseInt(document.getElementById('numBinsInput').value);

    if (!file) {
        document.getElementById('error').textContent = 'Please select a file first.';
        return;
    }

    if (isNaN(rMax) || rMax <= 0) {
        document.getElementById('error').textContent = 'Please enter a valid radius (rMax > 0).';
        return;
    }

    if (isNaN(numBins) || numBins < 1) {
        document.getElementById('error').textContent = 'Please enter a valid number of bins (numBins ≥ 1).';
        return;
    }

    readXYZFile(file, (positions, boxSize) => {
        const { theta, adf } = computeADF(positions, boxSize, rMax, numBins);
        plotADF(theta, adf);
        document.getElementById('error').textContent = '';
    });
});
