// --- Basic Setup ---
const audioContext = new (window.AudioContext || window.webkitAudioContext)();
let isAudioContextResumed = false;
let currentEngine = 'string'; // 'string', 'modal', 'hybrid'

// --- DOM Element References ---
const engineStringRadio = document.getElementById('engineString');
const engineModalRadio = document.getElementById('engineModal');
const engineHybridRadio = document.getElementById('engineHybrid');

// Modal / Hybrid Controls
const modalPresetASelect = document.getElementById('modalPresetA');
const modalPresetBSelect = document.getElementById('modalPresetB');
const modalMorphSlider = document.getElementById('modalMorph');
const chaosFeedbackSlider = document.getElementById('chaosFeedback');
const modalSpecificDiv = document.querySelector('.modal-hybrid-controls'); // Container div
const resonatorLegend = document.getElementById('resonator-legend'); // For dynamic title update


const dampingLabel = document.getElementById('damping-label');
const dampingHelp = document.getElementById('damping-help');
const pluckinessLabel = document.getElementById('pluckiness-label');
const pluckinessHelp = document.getElementById('pluckiness-help');

// Slider & General Control References
const frequencySlider = document.getElementById('frequency');
const dampingSlider = document.getElementById('damping');
const pluckinessSlider = document.getElementById('pluckiness');
const delayTimeSlider = document.getElementById('delayTime');
const delayFeedbackSlider = document.getElementById('delayFeedback');
const filterCutoffSlider = document.getElementById('filterCutoff');
const masterGainSlider = document.getElementById('masterGain');
const pluckButton = document.getElementById('pluckButton');
const controlsDiv = document.querySelector('.controls');

// NEW Effect Control References
const panLfoRateSlider = document.getElementById('panLfoRate');
const panLfoDepthSlider = document.getElementById('panLfoDepth');
const grainAmountSlider = document.getElementById('grainAmount');

// Value display spans
const freqValueSpan = document.getElementById('freq-value');
const dampingValueSpan = document.getElementById('damping-value');
const pluckinessValueSpan = document.getElementById('pluckiness-value');
const modalMorphValueSpan = document.getElementById('modalMorph-value'); // Already exists
const chaosFeedbackValueSpan = document.getElementById('chaosFeedback-value'); // Already exists
const delayTimeValueSpan = document.getElementById('delayTime-value');
const delayFeedbackValueSpan = document.getElementById('delayFeedback-value');
const filterCutoffValueSpan = document.getElementById('filterCutoff-value');
const masterGainValueSpan = document.getElementById('masterGain-value');
// NEW Value Spans
const panLfoRateValueSpan = document.getElementById('panLfoRate-value');
const panLfoDepthValueSpan = document.getElementById('panLfoDepth-value');
const grainAmountValueSpan = document.getElementById('grainAmount-value');

// NEW Visualizer Canvas References
const oscilloscopeCanvas = document.getElementById('oscilloscopeCanvas');
const spectrumCanvas = document.getElementById('spectrumCanvas');
let oscilloscopeCtx, spectrumCtx; // Will be assigned later


// Mouse Wheel Control Target
const mouseWheelControlTargetSelect = document.getElementById('mouseWheelControlTarget');


// --- Polyphony and State Management ---
const MAX_VOICES = 8; // Increased polyphony slightly
const activeVoices = []; // { key, type, nodes, startTime, pannerNode, outputNode }
const pressedKeys = new Set();

// --- Keyboard Mapping & Tuning --- (Standard)
const REFERENCE_KEY = 'k'; // A4
const keyToNoteMap = { /* ... (Standard mapping, unchanged) ... */
    'a': -12, 's': -10, 'd': -9, 'f': -7, 'g': -5, 'h': -4, 'j': -2, 'k': 0, 'l': 2, ';': 3,
    'w': -11, 'e': -8,  't': -6, 'y': -3, 'u': -1, 'o': 1,  'p': 4
};
function getFrequencyFromKey(key) {
    const keyLower = key.toLowerCase();
    const semitoneDifference = keyToNoteMap[keyLower];
    if (semitoneDifference === undefined) return null;
    const referenceFrequency = parseFloat(frequencySlider.value);
    return referenceFrequency * Math.pow(2, semitoneDifference / 12);
}


// --- Modal Synthesis Data --- (Added Levitation Chamber)
const modalPresets = {
    bell: { modes: [ { ratio: 0.56, gain: 0.4, decay: 1.8 }, { ratio: 0.92, gain: 0.25, decay: 2.5 }, { ratio: 1.19, gain: 0.15, decay: 2.2 }, { ratio: 1.71, gain: 0.1, decay: 1.5 }, { ratio: 2.00, gain: 0.3, decay: 1.9 }, { ratio: 2.74, gain: 0.08, decay: 1.2 }, { ratio: 3.00, gain: 0.07, decay: 1.4 }, { ratio: 3.76, gain: 0.05, decay: 0.9 }, ], decayScale: 4.0 },
    bar: { modes: [ { ratio: 1.00, gain: 0.6, decay: 2.0 }, { ratio: 2.76, gain: 0.2, decay: 1.5 }, { ratio: 5.40, gain: 0.1, decay: 1.0 }, { ratio: 8.93, gain: 0.05, decay: 0.7 }, { ratio: 13.34, gain: 0.03, decay: 0.5 }, ], decayScale: 3.0 },
    glass: { modes: [ { ratio: 1.00, gain: 0.5, decay: 1.5 }, { ratio: 1.61, gain: 0.3, decay: 1.8 }, { ratio: 2.41, gain: 0.15, decay: 1.2 }, { ratio: 3.88, gain: 0.1, decay: 0.9 }, { ratio: 5.11, gain: 0.08, decay: 0.7 }, { ratio: 7.23, gain: 0.05, decay: 0.5 }, ], decayScale: 2.5 },
    tube: { modes: [ { ratio: 1.00, gain: 0.6, decay: 1.8 }, { ratio: 2.00, gain: 0.3, decay: 1.5 }, { ratio: 3.00, gain: 0.15, decay: 1.2 }, { ratio: 4.00, gain: 0.1, decay: 0.9 }, { ratio: 5.00, gain: 0.05, decay: 0.6 }, ], decayScale: 3.5 },
    // --- NEW RESONATOR PRESET ---
    levitationChamber: {
        modes: [
            { ratio: 0.88, gain: 0.5, decay: 1.2 }, // Slightly flat fundamental, sharp decay
            { ratio: 1.00, gain: 0.10, decay: 4.0 }, // Persistent 'object' fundamental
            { ratio: 1.63, gain: 0.3, decay: 0.9 }, // Mode doublet around 1.65
            { ratio: 1.67, gain: 0.25, decay: 0.8 }, // Mode doublet partner
            { ratio: 2.91, gain: 0.15, decay: 1.1 }, // Higher non-harmonic partial
            { ratio: 4.05, gain: 0.08, decay: 0.6 }, // Even higher, fast decay
            { ratio: 5.52, gain: 0.05, decay: 0.4}   // Hint of another fast partial
        ],
        decayScale: 2.2 // Relatively short overall sustain compared to Bell/Tube
    }
};


// --- Audio Routing Setup --- (Added LFOs)
const summingBus = audioContext.createGain();
const masterVolume = audioContext.createGain();
const masterFilter = audioContext.createBiquadFilter();
const delay = audioContext.createDelay(1.0); // Max delay time
const delayFeedback = audioContext.createGain();
const delayWetGain = audioContext.createGain(); // Controls overall delay level
const delayInputGain = audioContext.createGain(); // Input level to delay line

// NEW LFOs and Controllers
const panLfo = audioContext.createOscillator();
const panLfoDepth = audioContext.createGain();
const delayLfo = audioContext.createOscillator();
const delayLfoDepth = audioContext.createGain(); // Controls grain modulation amount


// NEW Analyser Node
const analyserNode = audioContext.createAnalyser();
analyserNode.fftSize = 2048; // Power of 2. Determines resolution.
analyserNode.smoothingTimeConstant = 0.85; // 0 to 1. Smoothes frequency data.

// Data arrays for visualization
const analyserFrequencyCount = analyserNode.frequencyBinCount; // fftSize / 2
const timeDataArray = new Uint8Array(analyserFrequencyCount); // For oscilloscope
const frequencyDataArray = new Uint8Array(analyserFrequencyCount); // For spectrum

// --- Connections (REVISED) ---
summingBus.connect(masterFilter); // Direct signal to filter
summingBus.connect(delayInputGain); // Direct signal to delay input

// Delay Path
delayInputGain.connect(delay);
delay.connect(delayFeedback);
delayFeedback.connect(delay);
delay.connect(delayWetGain);
delayWetGain.connect(masterFilter); // Delayed signal ALSO goes to filter

// Filter connects to Analyser
masterFilter.connect(analyserNode);

// Analyser connects to Master Volume
analyserNode.connect(masterVolume);

// Master Volume connects to Output
masterVolume.connect(audioContext.destination);

// LFO Connections (remain the same, not part of main audio path)
panLfo.connect(panLfoDepth);
// panLfoDepth connects dynamically to each voice's pannerNode.pan
delayLfo.connect(delayLfoDepth);
delayLfoDepth.connect(delay.delayTime);


// Initial Effect Parameter Settings
masterFilter.type = 'lowpass';
masterFilter.Q.value = 1;
delayWetGain.gain.value = 0.7; // Overall delay mix
delayInputGain.gain.value = 0.6; // Signal fed into delay

// Start LFOs
panLfo.type = 'sine';
panLfo.frequency.value = parseFloat(panLfoRateSlider.value); // Initial rate
panLfoDepth.gain.value = parseFloat(panLfoDepthSlider.value); // Initial depth (0 = no panning)
panLfo.start();

delayLfo.type = 'triangle'; // Triangle or Sine often work well for smooth delay modulation
delayLfo.frequency.value = 6; // A typical rate for granular feel (adjust to taste)
delayLfoDepth.gain.value = parseFloat(grainAmountSlider.value); // Initial grain amount (0 = no mod)
delayLfo.start();



// --- Synthesis Functions ---

// Plucked String (Karplus-Strong) - Unchanged
function createPluckedStringVoice(frequency, damping, pluckiness) {
    if (!audioContext || audioContext.state !== 'running') return null;
    const sampleRate = audioContext.sampleRate;
    const N = Math.round(sampleRate / frequency);
    const duration = 3.5; // Approx max duration
    const bufferSize = Math.max(N, sampleRate * duration); // Ensure buffer is large enough
    const buffer = audioContext.createBuffer(1, bufferSize, sampleRate);
    const data = buffer.getChannelData(0);

    // Create noise burst
    const delayLine = new Float32Array(N);
    for (let i = 0; i < N; i++) {
        delayLine[i] = Math.random() * 2 - 1;
    }

    // Simple low-pass filter on noise based on pluckiness
    const pluckFilterFactor = 1.0 - pluckiness; // How much of previous sample to mix in
    let prevNoise = 0;
    for (let i = 0; i < N; i++) {
        let currentSample = delayLine[i];
        // Simple 1-pole LPF: y[n] = x[n]*pluckiness + y[n-1]*(1-pluckiness)
        delayLine[i] = currentSample * pluckiness + prevNoise * pluckFilterFactor;
        prevNoise = delayLine[i];
    }

    // Karplus-Strong algorithm
    let currentIndex = 0;
    for (let i = 0; i < bufferSize; i++) {
        const currentSample = delayLine[currentIndex % N];
        const nextIndex = (currentIndex + 1) % N;
        const nextSample = delayLine[nextIndex];
        // Averaging filter + damping
        const newSample = (currentSample + nextSample) * 0.5 * damping;
        data[i] = currentSample; // Write current sample to output buffer
        delayLine[currentIndex % N] = newSample; // Update delay line
        currentIndex++;
    }

    const sourceNode = audioContext.createBufferSource();
    sourceNode.buffer = buffer;

    const gainNode = audioContext.createGain();
    gainNode.gain.setValueAtTime(0.8, audioContext.currentTime); // Initial gain

    sourceNode.connect(gainNode);

    return { type: 'string', outputNode: gainNode, stopNodes: [sourceNode] };
}


// Modal Synthesis - Morphing & Chaos - Unchanged from previous combined code
function createModalVoice(frequency, presetAData, presetBData, morphValue, globalDampingValue, chaosFeedbackAmount) {
    if (!audioContext || audioContext.state !== 'running') return null;
    if (!presetAData || !presetBData) { console.warn("Missing preset data for modal voice"); return null; }

    const sampleRate = audioContext.sampleRate;
    const now = audioContext.currentTime;
    const voiceOutputGain = audioContext.createGain();
    voiceOutputGain.gain.setValueAtTime(1.0, now);

    const oscillators = [];
    const gainNodes = []; // Store primary gain nodes for feedback access
    const feedbackGainNodes = [];
    let longestDecay = 0;

    // Interpolate decay scale and calculate multiplier based on global damping
    const interpDecayScale = presetAData.decayScale + (presetBData.decayScale - presetAData.decayScale) * morphValue;
    // Inverse relationship for modal: Lower damping slider value = longer decay
    const decayMultiplier = Math.max(0.1, interpDecayScale * (1.05 - globalDampingValue));

    const numModes = Math.max(presetAData.modes.length, presetBData.modes.length);

    for (let i = 0; i < numModes; i++) {
        const modeA = presetAData.modes[i];
        const modeB = presetBData.modes[i];

        let interpRatio, interpGain, interpDecay, effectiveGain;

        // Interpolate mode parameters (ratio, gain, decay)
        if (modeA && modeB) {
            interpRatio = modeA.ratio + (modeB.ratio - modeA.ratio) * morphValue;
            interpGain = modeA.gain + (modeB.gain - modeA.gain) * morphValue;
            interpDecay = modeA.decay + (modeB.decay - modeA.decay) * morphValue;
            effectiveGain = interpGain;
        } else if (modeA) { // Only exists in A, fade out
            interpRatio = modeA.ratio; interpGain = modeA.gain; interpDecay = modeA.decay;
            effectiveGain = interpGain * (1.0 - morphValue);
        } else if (modeB) { // Only exists in B, fade in
            interpRatio = modeB.ratio; interpGain = modeB.gain; interpDecay = modeB.decay;
            effectiveGain = interpGain * morphValue;
        } else {
            continue; // Skip if mode doesn't exist in either
        }

        const modeFreq = frequency * interpRatio;
        const decayTime = interpDecay * decayMultiplier; // Apply calculated decay

        // Skip invalid or inaudible modes
        if (modeFreq > sampleRate / 2 || modeFreq <= 0 || effectiveGain < 0.001 || decayTime <= 0.001) continue;

        const osc = audioContext.createOscillator();
        const gain = audioContext.createGain(); // This mode's primary gain envelope

        osc.frequency.setValueAtTime(modeFreq, now);
        osc.type = 'sine';

        // Gain envelope: quick attack, exponential decay
        gain.gain.setValueAtTime(0, now);
        gain.gain.linearRampToValueAtTime(effectiveGain, now + 0.005); // Short attack
        gain.gain.setTargetAtTime(0.0001, now + 0.01, decayTime / 4); // Exponential decay target

        osc.connect(gain);
        gain.connect(voiceOutputGain);

        osc.start(now);
        const stopTime = now + decayTime * 1.5 + 0.2; // Schedule stop well after decay finishes
        osc.stop(stopTime);
        if (decayTime > longestDecay) longestDecay = decayTime;

        oscillators.push(osc);
        gainNodes.push(gain); // Store for chaos feedback
    }

    // --- Chaotic Feedback Implementation ---
    const CHAOS_FEEDBACK_SCALING = 0.15; // Safety scaling factor
    if (chaosFeedbackAmount > 0 && gainNodes.length > 1) {
        gainNodes.forEach((sourceGainNode, index) => {
            // Feed back to the *next* mode's gain PARAMETER (wraps around)
            const targetGainNode = gainNodes[(index + 1) % gainNodes.length];
            const feedbackGain = audioContext.createGain();
            const feedbackLevel = chaosFeedbackAmount * CHAOS_FEEDBACK_SCALING;
            feedbackGain.gain.setValueAtTime(feedbackLevel, now);

            sourceGainNode.connect(feedbackGain); // Connect source output
            feedbackGain.connect(targetGainNode.gain); // Connect TO the gain AudioParam

            feedbackGainNodes.push(feedbackGain); // Keep track if needed later
        });
        // console.log(`Chaos Feedback enabled for modal voice. Level: ${feedbackLevel.toFixed(3)}`);
    }

    if (oscillators.length === 0) { console.warn("No valid modes created for modal voice"); return null;}

    return {
        type: 'modal',
        outputNode: voiceOutputGain,
        stopNodes: oscillators,
        feedbackNodes: feedbackGainNodes // Include if cleanup is ever needed
    };
}

// Hybrid Resonance - Morphing, Chaos & Virtual String Excit. - Unchanged from prev.
function createHybridResonanceVoice(frequency, presetAData, presetBData, morphValue, excitationPluckiness, modalDampingValue, chaosFeedbackAmount) {
     if (!audioContext || audioContext.state !== 'running') return null;
     if (!presetAData || !presetBData) { console.warn("Missing preset data for hybrid voice"); return null; }

     const sampleRate = audioContext.sampleRate;
     const now = audioContext.currentTime;
     const voiceOutputGain = audioContext.createGain();
     voiceOutputGain.gain.setValueAtTime(1.0, now); // Start full volume, envelopes shape sound

     const oscillators = [];
     const primaryGainNodes = []; // Store base decay gain nodes for feedback
     const feedbackGainNodes = [];
     let longestDecay = 0;

     // --- Virtual String Excitation Simulation ---
     const virtualStringFundamental = frequency;
     const VIRTUAL_STRING_BASE_DECAY = 1.5; // Base decay of fundamental (seconds)
     const VIRTUAL_STRING_DAMPING_PER_OCTAVE = 0.5; // Higher harmonics decay faster
     const MAX_VIRTUAL_HARMONICS = 20;
     const virtualHarmonicAmplitudes = []; // Calculated amplitude for each harmonic
     const baseAmplitudes = []; let amplitudeSum = 0;

     // Initial amplitudes (1/n)
     for(let i = 1; i <= MAX_VIRTUAL_HARMONICS; i++) { baseAmplitudes[i] = 1 / i; }

     // Simulate pluckiness filtering on harmonic amplitudes (simple LPF)
     const filterFactor = 1.0 - excitationPluckiness; let previousFilteredAmp = 0;
     for(let i = 1; i <= MAX_VIRTUAL_HARMONICS; i++) {
         const filteredAmp = baseAmplitudes[i] * excitationPluckiness + previousFilteredAmp * filterFactor;
         virtualHarmonicAmplitudes[i] = Math.max(0, Math.min(1.0, filteredAmp));
         previousFilteredAmp = filteredAmp;
         amplitudeSum += virtualHarmonicAmplitudes[i];
     }
     // Normalize amplitudes roughly (adjust targetSum to control overall brightness)
     if (amplitudeSum > 0) {
         const targetSum = 2.0; // Controls overall brightness/level from excitation
         const normFactor = targetSum / amplitudeSum;
         for(let i = 1; i <= MAX_VIRTUAL_HARMONICS; i++) {
             virtualHarmonicAmplitudes[i] = (virtualHarmonicAmplitudes[i] || 0) * normFactor;
         }
     } else {
         for(let i = 1; i <= MAX_VIRTUAL_HARMONICS; i++) { virtualHarmonicAmplitudes[i] = 0;}
     }
     // --- End Virtual String ---


     // --- Modal Resonator Setup with Morphing ---
     const interpDecayScale = presetAData.decayScale + (presetBData.decayScale - presetAData.decayScale) * morphValue;
     const decayMultiplier = Math.max(0.1, interpDecayScale * (1.05 - modalDampingValue)); // Inverse relationship again

     const numModes = Math.max(presetAData.modes.length, presetBData.modes.length);

     for (let i = 0; i < numModes; i++) {
         const modeA = presetAData.modes[i];
         const modeB = presetBData.modes[i];
         let interpRatio, interpGain, interpDecay, effectiveBaseGain;

         // Interpolate Mode Parameters (Same logic as in createModalVoice)
         if (modeA && modeB) {
              interpRatio = modeA.ratio + (modeB.ratio - modeA.ratio) * morphValue;
              interpGain = modeA.gain + (modeB.gain - modeA.gain) * morphValue;
              interpDecay = modeA.decay + (modeB.decay - modeA.decay) * morphValue;
              effectiveBaseGain = interpGain;
          } else if (modeA) {
              interpRatio = modeA.ratio; interpGain = modeA.gain; interpDecay = modeA.decay;
              effectiveBaseGain = interpGain * (1.0 - morphValue);
          } else if (modeB) {
              interpRatio = modeB.ratio; interpGain = modeB.gain; interpDecay = modeB.decay;
              effectiveBaseGain = interpGain * morphValue;
          } else {
              continue;
          }
          // --- End Interpolation ---

         const osc = audioContext.createOscillator();
         const primaryModeGain = audioContext.createGain(); // Controls the modal decay envelope
         const activationGain = audioContext.createGain(); // Controls amplitude based on string excitation

         const modeFreq = frequency * interpRatio;
         const baseModeDecayTime = interpDecay * decayMultiplier; // Base decay for this mode

         if (modeFreq > sampleRate / 2 || modeFreq <= 0 || effectiveBaseGain < 0.001 || baseModeDecayTime <= 0.001) continue;

         osc.frequency.setValueAtTime(modeFreq, now);
         osc.type = 'sine';

         // --- Calculate Virtual String Influence ---
         let closestHarmonicIdx = -1; let minFreqDiff = Infinity;
         // Find the virtual string harmonic closest to this mode's frequency
         for (let h = 1; h <= MAX_VIRTUAL_HARMONICS; h++) {
             const harmonicFreq = virtualStringFundamental * h;
             const freqDiff = Math.abs(modeFreq - harmonicFreq);
             if (freqDiff < minFreqDiff) { minFreqDiff = freqDiff; closestHarmonicIdx = h; }
         }
         // Calculate activation based on frequency proximity and harmonic amplitude
         let harmonicActivation = 0;
         if (closestHarmonicIdx > 0) {
             // Proximity factor based on semitone difference (tune multiplier '2.0' for sensitivity)
             const semitoneDiff = 12 * Math.log2(modeFreq / (virtualStringFundamental * closestHarmonicIdx));
             const proximityFactor = Math.max(0, 1.0 - Math.abs(semitoneDiff) * 2.0);
             harmonicActivation = proximityFactor * (virtualHarmonicAmplitudes[closestHarmonicIdx] || 0);
         }
         // --- End String Influence ---

         // --- Set up Gain Nodes using *interpolated* base gain & activation ---
         // Primary Gain: Base modal decay, scaled slightly by activation
         primaryModeGain.gain.setValueAtTime(0, now);
         // Initial gain includes base interpolated gain, boosted by activation
         primaryModeGain.gain.linearRampToValueAtTime(effectiveBaseGain * Math.max(0.1, harmonicActivation + 0.1), now + 0.010);
         primaryModeGain.gain.setTargetAtTime(0.0001, now + 0.01, baseModeDecayTime / 4); // Modal decay

         // Activation Gain: Faster decay based on the virtual harmonic's decay
         const harmonicDecayTime = VIRTUAL_STRING_BASE_DECAY * Math.pow(VIRTUAL_STRING_DAMPING_PER_OCTAVE, Math.log2(Math.max(1, closestHarmonicIdx)));
         activationGain.gain.setValueAtTime(Math.max(0.0001, harmonicActivation), now); // Set initial level based on activation
         activationGain.gain.setTargetAtTime(0.0001, now + 0.005, Math.max(0.01, harmonicDecayTime / 3)); // Faster decay

         // Connect nodes: osc -> primaryGain (modal decay) -> activationGain (string impulse) -> voiceOutput
         osc.connect(primaryModeGain);
         primaryModeGain.connect(activationGain);
         activationGain.connect(voiceOutputGain);

         // --- Start & Schedule Stop ---
         osc.start(now);
         // Stop time should accommodate the longer of the modal or harmonic decay
         const stopTime = now + Math.max(baseModeDecayTime * 1.5, harmonicDecayTime * 1.5) + 0.2;
         osc.stop(stopTime);
         if (baseModeDecayTime > longestDecay) longestDecay = baseModeDecayTime;

         oscillators.push(osc);
         primaryGainNodes.push(primaryModeGain); // Store for potential feedback
     }

    // --- Chaotic Feedback Implementation --- (Connects primaryGain nodes)
    const CHAOS_FEEDBACK_SCALING = 0.15;
    if (chaosFeedbackAmount > 0 && primaryGainNodes.length > 1) {
          primaryGainNodes.forEach((sourceGainNode, index) => {
             const targetGainNode = primaryGainNodes[(index + 1) % primaryGainNodes.length];
             const feedbackGain = audioContext.createGain();
             const feedbackLevel = chaosFeedbackAmount * CHAOS_FEEDBACK_SCALING;
             feedbackGain.gain.setValueAtTime(feedbackLevel, now);

             sourceGainNode.connect(feedbackGain); // From source output
             feedbackGain.connect(targetGainNode.gain); // TO target gain param

             feedbackGainNodes.push(feedbackGain);
          });
        //   console.log(`Chaos Feedback enabled for hybrid voice. Level: ${feedbackLevel.toFixed(3)}`);
    }


    if (oscillators.length === 0) { console.warn(`No valid modes/harmonics for hybrid. Freq: ${frequency.toFixed(1)}, Morph:${morphValue.toFixed(2)}`); return null; }

    return {
        type: 'hybrid',
        outputNode: voiceOutputGain,
        stopNodes: oscillators,
        feedbackNodes: feedbackGainNodes
    };
}

// --- Function to Start a Voice --- MODIFIED FOR LFO PANNING
function playNote(key, frequency) {
     if (!isAudioContextResumed || audioContext.state !== 'running') {
        console.log("Audio Context not running."); if (audioContext.state === 'suspended') { resumeAudioContext().then(() => { if (audioContext.state === 'running') playNote(key, frequency); }); } return;
     }

    // --- Voice Management (Stealing) ---
    if (activeVoices.length >= MAX_VOICES) {
        activeVoices.sort((a, b) => a.startTime - b.startTime); // Find oldest voice
        const stolenVoice = activeVoices.shift(); // Remove it
        console.log(`Max voices reached. Stealing oldest voice (Key: ${stolenVoice.key}, Type: ${stolenVoice.type}).`);
        // Fade out and stop the stolen voice quickly
        stolenVoice.outputNode.gain.cancelScheduledValues(audioContext.currentTime);
        stolenVoice.outputNode.gain.setTargetAtTime(0.0, audioContext.currentTime, 0.015);
        stolenVoice.nodes.stopNodes.forEach(node => { try { node.stop(audioContext.currentTime + 0.05); } catch(e) {/* Might already be stopped */} });
        try { // Disconnect explicitly to be sure
            panLfoDepth.disconnect(stolenVoice.pannerNode.pan); // Disconnect LFO from this specific panner param
            stolenVoice.pannerNode.disconnect();
        } catch(e) { /* Already disconnected maybe */ }
        pressedKeys.delete(stolenVoice.key); // Allow retriggering the key
    }

    // --- Create Voice Based on Engine ---
    let voiceData = null;
    const now = audioContext.currentTime;

    if (currentEngine === 'string') {
        const damping = parseFloat(dampingSlider.value);
        const pluckiness = parseFloat(pluckinessSlider.value);
        voiceData = createPluckedStringVoice(frequency, damping, pluckiness);
    } else { // Modal or Hybrid
        const presetAName = modalPresetASelect.value;
        const presetBName = modalPresetBSelect.value;
        const presetAData = modalPresets[presetAName];
        const presetBData = modalPresets[presetBName];
        const morphValue = parseFloat(modalMorphSlider.value);
        const dampingValue = parseFloat(dampingSlider.value); // Read as "Resonance" for Modal/Hybrid
        const chaosValue = parseFloat(chaosFeedbackSlider.value);

        if (!presetAData || !presetBData) {
             console.error("Could not find preset data for morphing!");
             return; // Stop if presets are missing
        }

        if (currentEngine === 'modal') {
            voiceData = createModalVoice(frequency, presetAData, presetBData, morphValue, dampingValue, chaosValue);
        } else { // Hybrid engine
            const pluckinessValue = parseFloat(pluckinessSlider.value); // Read as "Excitation Brightness"
            voiceData = createHybridResonanceVoice(frequency, presetAData, presetBData, morphValue, pluckinessValue, dampingValue, chaosValue);
        }
    }


    // --- Connect and Store Voice ---
    if (voiceData) {
        const pannerNode = audioContext.createStereoPanner();
        // REMOVED Random Panning: pannerNode.pan.setValueAtTime((Math.random() * 1.6) - 0.8, now);
        // Connect the GLOBAL Pan LFO Depth controller to THIS voice's panner parameter
        panLfoDepth.connect(pannerNode.pan);

        voiceData.outputNode.connect(pannerNode);
        pannerNode.connect(summingBus); // Voice Output -> Panner -> Summing Bus

        if (voiceData.type === 'string' && voiceData.stopNodes && voiceData.stopNodes[0]) {
             voiceData.stopNodes[0].start(now); // Start buffer source immediately
             // Cleanup for buffer sources needs onended
             voiceData.stopNodes[0].onended = () => cleanupVoice(voiceData.stopNodes[0]);
        }
        // Modal & Hybrid oscillators rely on scheduled stop() or voice stealing for cleanup

         activeVoices.push({
             key: key,
             type: voiceData.type,
             nodes: voiceData,
             outputNode: voiceData.outputNode, // Main gain node for the voice
             pannerNode: pannerNode,           // Panner specific to this voice
             startTime: now
         });
     } else {
         console.log("Failed to create voice data for key:", key);
     }
}

// Helper function to remove voice from active list (mostly for string buffer sources)
function cleanupVoice(stoppedNode) {
    // Find the voice associated with the stopped node
    const voiceIndex = activeVoices.findIndex(v => v.nodes && v.nodes.stopNodes && v.nodes.stopNodes.includes(stoppedNode));
    if (voiceIndex > -1) {
        const removedVoice = activeVoices.splice(voiceIndex, 1)[0];
        // Disconnect nodes to release resources, although garbage collection should handle it eventually
        try { removedVoice.outputNode.disconnect(); } catch(e) {}
        try {
            panLfoDepth.disconnect(removedVoice.pannerNode.pan); // Disconnect LFO from this panner
            removedVoice.pannerNode.disconnect();
        } catch(e) {}
        // console.log(`Cleaned up voice for key: ${removedVoice.key}`);
    }
}


// Function to Resume Audio Context (Standard)
async function resumeAudioContext() {
    if (audioContext.state === 'suspended') {
        console.log('Attempting to resume AudioContext...');
        await audioContext.resume();
    }
    if (audioContext.state === 'running') {
        isAudioContextResumed = true;
        console.log('AudioContext is running.');
    } else {
        console.warn('AudioContext could not be resumed. State:', audioContext.state);
    }
}

// Update UI based on Engine - Handles dynamic labels/ranges/visibility
function updateUIForEngine() {
    const body = document.body;

    // --- Step 1: Clear existing engine classes ---
    body.classList.remove('engine-string', 'engine-modal', 'engine-hybrid');

    // --- Step 2: Get references needed for conditional logic ---
    const modalHybridWheelOptions = document.querySelectorAll('#mouseWheelControlTarget .modal-hybrid-wheel');

    // --- Step 3: Apply engine-specific settings ---
    if (currentEngine === 'string') {
        body.classList.add('engine-string');

        // Update Labels & Help Text
        dampingLabel.textContent = "Damping (Decay):";
        dampingHelp.textContent = "(Higher = longer decay)";
        pluckinessLabel.textContent = "Pluck Stiffness:";
        pluckinessHelp.textContent = "(Higher = brighter)";
        if (resonatorLegend) resonatorLegend.textContent = "Resonator Controls"; // Reset (might be hidden)

        // Update Damping Slider Range & Value (for Karplus-Strong)
        dampingSlider.min = "0.95";
        dampingSlider.max = "0.999";
        dampingSlider.step = "0.001";
        // Coerce value into new range if needed
        if (parseFloat(dampingSlider.value) < 0.95 || parseFloat(dampingSlider.value) > 0.999 || isNaN(parseFloat(dampingSlider.value))) {
            dampingSlider.value = 0.996; // Default string damping
        }

        // Update Mouse Wheel Options Visibility
        modalHybridWheelOptions.forEach(opt => opt.style.display = 'none');
        // Switch wheel target if current selection is now hidden
        if (['modalMorph', 'chaosFeedback'].includes(mouseWheelControlTargetSelect.value)) {
            mouseWheelControlTargetSelect.value = 'frequency'; // Default safe choice
        }

    } else { // Modal or Hybrid (Share many UI properties)
        // Set body class based on specific engine
        body.classList.add(currentEngine === 'modal' ? 'engine-modal' : 'engine-hybrid');

        // Update Labels & Help Text
        dampingLabel.textContent = "Resonance Decay:";
        dampingHelp.textContent = "(Lower = longer)";
        if (resonatorLegend) resonatorLegend.textContent = currentEngine === 'modal' ? "Modal Resonator" : "Hybrid Resonator";

        if (currentEngine === 'hybrid') {
            pluckinessLabel.textContent = "Excitation Brightness:"; // Label specific to Hybrid
            pluckinessHelp.textContent = "(Virtual string excitation)";
        } else { // Modal - Pluckiness control is hidden via CSS
             pluckinessLabel.textContent = "Pluck/Excitation:"; // Reset label (it's hidden)
             pluckinessHelp.textContent = "";
        }


        // Update Damping Slider Range & Value (for Modal/Hybrid resonance)
        dampingSlider.min = "0.1";
        dampingSlider.max = "1.0";
        dampingSlider.step = "0.01";
        // Coerce value into new range if needed
        if (parseFloat(dampingSlider.value) > 1.0 || parseFloat(dampingSlider.value) < 0.1 || isNaN(parseFloat(dampingSlider.value))) {
            dampingSlider.value = 0.5; // Default mid-value for modal/hybrid
        }

        // Update Mouse Wheel Options Visibility
         modalHybridWheelOptions.forEach(opt => opt.style.display = 'block'); // Show morph/chaos
    }

    // --- Step 4: Update value displays AFTER potentially changing slider values ---
    // Format damping value based on range/step
    const dampVal = parseFloat(dampingSlider.value);
    dampingValueSpan.textContent = (currentEngine === 'string' ? dampVal.toFixed(3) : dampVal.toFixed(2));

    // Update pluckiness display only if its control is relevant for the current engine
    if (currentEngine === 'string' || currentEngine === 'hybrid') {
        pluckinessValueSpan.textContent = parseFloat(pluckinessSlider.value).toFixed(2);
    } else {
        pluckinessValueSpan.textContent = "---"; // Indicate N/A when hidden
    }

    // --- Step 5: Rely on CSS rules linked to body classes ---
    //     (e.g., body.engine-string .string-specific { display: flex; })
    //     to handle showing/hiding fieldsets/control groups.
}


// --- Event Listeners ---

// Engine Switch
engineStringRadio.addEventListener('change', () => { if (engineStringRadio.checked) { currentEngine = 'string'; updateUIForEngine(); } });
engineModalRadio.addEventListener('change', () => { if (engineModalRadio.checked) { currentEngine = 'modal'; updateUIForEngine(); } });
engineHybridRadio.addEventListener('change', () => { if (engineHybridRadio.checked) { currentEngine = 'hybrid'; updateUIForEngine(); } });


// Pluck Button
pluckButton.addEventListener('click', async () => {
    await resumeAudioContext(); // Ensure context is running
    if (isAudioContextResumed) {
        const referenceFrequency = parseFloat(frequencySlider.value);
        playNote(REFERENCE_KEY, referenceFrequency); // Trigger note using reference key
    }
});

// Keyboard
window.addEventListener('keydown', async (event) => {
    // Ignore if modifier keys are pressed or if it's a repeat event
    if (event.metaKey || event.ctrlKey || event.altKey || event.repeat) return;
    const frequency = getFrequencyFromKey(event.key);
    if (frequency !== null) { // Check if the key is mapped
        event.preventDefault(); // Prevent default actions like typing
        if (audioContext.state === 'suspended') await resumeAudioContext();
        const keyLower = event.key.toLowerCase();
        // Only play if the key isn't already pressed (prevents retriggering noise)
        if (!pressedKeys.has(keyLower)) {
            pressedKeys.add(keyLower);
            playNote(keyLower, frequency);
        }
    }
});
window.addEventListener('keyup', (event) => {
    const keyLower = event.key.toLowerCase();
    if (keyToNoteMap.hasOwnProperty(keyLower)) {
        pressedKeys.delete(keyLower); // Remove key from pressed set
        // Note: We don't explicitly stop the sound here; relying on decay or voice stealing.
        // Could implement note-off messages later if needed.
    }
});


// --- Slider Updates --- (Connect UI to Audio Params)
frequencySlider.oninput = () => freqValueSpan.textContent = frequencySlider.value;
dampingSlider.oninput = () => { // Display formatted value
    const dampVal = parseFloat(dampingSlider.value);
    dampingValueSpan.textContent = (currentEngine === 'string' ? dampVal.toFixed(3) : dampVal.toFixed(2));
};
pluckinessSlider.oninput = () => pluckinessValueSpan.textContent = parseFloat(pluckinessSlider.value).toFixed(2);
masterGainSlider.oninput = () => { masterGainValueSpan.textContent = masterGainSlider.value; masterVolume.gain.setTargetAtTime(parseFloat(masterGainSlider.value), audioContext.currentTime, 0.01); };
delayTimeSlider.oninput = () => { delayTimeValueSpan.textContent = delayTimeSlider.value; delay.delayTime.setTargetAtTime(parseFloat(delayTimeSlider.value), audioContext.currentTime, 0.01); };
delayFeedbackSlider.oninput = () => { delayFeedbackValueSpan.textContent = delayFeedbackSlider.value; delayFeedback.gain.setTargetAtTime(parseFloat(delayFeedbackSlider.value), audioContext.currentTime, 0.01); };
filterCutoffSlider.oninput = () => { filterCutoffValueSpan.textContent = filterCutoffSlider.value; masterFilter.frequency.setTargetAtTime(parseFloat(filterCutoffSlider.value), audioContext.currentTime, 0.01); };
modalMorphSlider.oninput = () => modalMorphValueSpan.textContent = parseFloat(modalMorphSlider.value).toFixed(2);
chaosFeedbackSlider.oninput = () => chaosFeedbackValueSpan.textContent = parseFloat(chaosFeedbackSlider.value).toFixed(2);

// NEW Slider Listeners
panLfoRateSlider.oninput = () => { panLfoRateValueSpan.textContent = parseFloat(panLfoRateSlider.value).toFixed(2); panLfo.frequency.setTargetAtTime(parseFloat(panLfoRateSlider.value), audioContext.currentTime, 0.01); };
panLfoDepthSlider.oninput = () => { panLfoDepthValueSpan.textContent = parseFloat(panLfoDepthSlider.value).toFixed(2); panLfoDepth.gain.setTargetAtTime(parseFloat(panLfoDepthSlider.value), audioContext.currentTime, 0.01); };
grainAmountSlider.oninput = () => { grainAmountValueSpan.textContent = parseFloat(grainAmountSlider.value).toFixed(3); delayLfoDepth.gain.setTargetAtTime(parseFloat(grainAmountSlider.value), audioContext.currentTime, 0.01); };


// --- Mouse Wheel Control ---
function handleWheelControl(event) {
    const targetId = mouseWheelControlTargetSelect.value;
    const slider = document.getElementById(targetId);

    if (!slider || slider.type !== 'range') {
        // Only proceed if the target is a valid range slider
        // Check if target is visible (don't control hidden sliders)
        const style = window.getComputedStyle(slider);
        if (!style || style.display === 'none' || slider.closest('.engine-specific:not([style*="display: flex"]):not([style*="display: block"])')) {
             // Heuristic check: if slider or its engine-specific parent is hidden
             // console.warn("Wheel target is hidden:", targetId);
            return;
        }
        //console.warn("Wheel target not found or not a range slider:", targetId);
        //return; // Commented out to allow page scroll if target isn't slider
    } else {
         // Prevent page scrolling ONLY when wheeling over a valid, visible slider target
         event.preventDefault();
         event.stopPropagation();
    }


    // Calculation Logic (Only if slider found and is range)
    if(slider && slider.type === 'range') {
        const delta = event.deltaY; // Positive for down/backward, Negative for up/forward
        const direction = delta < 0 ? 1 : -1; // 1 for increase, -1 for decrease

        const step = parseFloat(slider.step) || 0.01;
        const min = parseFloat(slider.min);
        const max = parseFloat(slider.max);
        const currentValue = parseFloat(slider.value);
        const range = max - min;

        // Adjust step size based on parameter sensitivity or slider range
        let effectiveStep = step;
        // Use larger steps for wide ranges like frequency or cutoff
        if (range > 1000) { effectiveStep = Math.max(step, range / 100); }
        else if (step < 0.01) { effectiveStep *= 3; } // More sensitivity for fine controls
        else if (step < 0.1) { effectiveStep *= 2; }


        let newValue = currentValue + direction * effectiveStep;

        // Clamp the value within min/max
        newValue = Math.max(min, Math.min(max, newValue));

        // Round to slider's step precision to avoid floating point issues
        const numDecimalPlaces = (step.toString().split('.')[1] || '').length;
        if (numDecimalPlaces > 0) {
            newValue = parseFloat(newValue.toFixed(numDecimalPlaces));
        } else {
            newValue = Math.round(newValue); // For integer steps
        }


        // Apply the change
        slider.value = newValue;

        // Trigger the 'input' event to update UI spans and audio params
        slider.dispatchEvent(new Event('input', { bubbles: true }));
    }
}

// Attach wheel listener to the main controls container
controlsDiv.addEventListener('wheel', handleWheelControl, { passive: false }); // passive: false needed for preventDefault


function drawOscilloscope() {
    if (!oscilloscopeCtx) return; // Skip if context not ready

    analyserNode.getByteTimeDomainData(timeDataArray); // Fill array with waveform data

    oscilloscopeCtx.fillStyle = 'rgb(44, 62, 80)'; // Match --input-bg approx
    oscilloscopeCtx.fillRect(0, 0, oscilloscopeCanvas.width, oscilloscopeCanvas.height);

    oscilloscopeCtx.lineWidth = 2;
    oscilloscopeCtx.strokeStyle = 'rgb(52, 152, 219)'; // Match --accent-color-primary

    oscilloscopeCtx.beginPath();

    const sliceWidth = oscilloscopeCanvas.width * 1.0 / analyserFrequencyCount;
    let x = 0;

    for (let i = 0; i < analyserFrequencyCount; i++) {
        const v = timeDataArray[i] / 128.0; // Normalize to 0.0 - 2.0
        const y = v * oscilloscopeCanvas.height / 2;

        if (i === 0) {
            oscilloscopeCtx.moveTo(x, y);
        } else {
            oscilloscopeCtx.lineTo(x, y);
        }

        x += sliceWidth;
    }

    oscilloscopeCtx.lineTo(oscilloscopeCanvas.width, oscilloscopeCanvas.height / 2);
    oscilloscopeCtx.stroke();
}

function drawSpectrum() {
    if (!spectrumCtx) return; // Skip if context not ready

    analyserNode.getByteFrequencyData(frequencyDataArray); // Fill array with frequency data

    spectrumCtx.fillStyle = 'rgb(44, 62, 80)'; // Match --input-bg approx
    spectrumCtx.fillRect(0, 0, spectrumCanvas.width, spectrumCanvas.height);

    const barWidth = (spectrumCanvas.width / analyserFrequencyCount) * 2.5; // Adjust multiplier for spacing
    let barHeight;
    let x = 0;

    for (let i = 0; i < analyserFrequencyCount; i++) {
        barHeight = frequencyDataArray[i]; // Value from 0-255

        // Use a gradient or accent color
        // spectrumCtx.fillStyle = 'rgb(' + (barHeight+100) + ',50,50)'; // Simple heatmap example
        spectrumCtx.fillStyle = 'rgb(46, 204, 113)'; // --accent-color-positive (greenish)

        // Map bar height (0-255) to canvas height
        const y = spectrumCanvas.height - (barHeight / 255 * spectrumCanvas.height);
        spectrumCtx.fillRect(x, y, barWidth, spectrumCanvas.height - y);

        x += barWidth + 1; // Add 1 for spacing between bars
    }
}

// --- Animation Loop ---
function visualize() {
    // Call the drawing functions
    drawOscilloscope();
    drawSpectrum();

    // Schedule the next frame
    requestAnimationFrame(visualize);
}


// --- Initial Setup ---
function initializeControlsAndAudio() {
    // Get canvas contexts
    if (oscilloscopeCanvas) oscilloscopeCtx = oscilloscopeCanvas.getContext('2d');
    if (spectrumCanvas) spectrumCtx = spectrumCanvas.getContext('2d');

     // Fix canvas HD-DPI scaling issues if contexts are available
     function setupCanvas(canvas, ctx) {
        if (!canvas || !ctx) return;
        const dpr = window.devicePixelRatio || 1;
        const rect = canvas.getBoundingClientRect();
        canvas.width = rect.width * dpr;
        canvas.height = rect.height * dpr;
        ctx.scale(dpr, dpr); // Scale drawings
        // Set width/height styles to ensure layout size is correct
        canvas.style.width = `${rect.width}px`;
        canvas.style.height = `${rect.height}px`;
    }
    setupCanvas(oscilloscopeCanvas, oscilloscopeCtx);
    setupCanvas(spectrumCanvas, spectrumCtx);


    // Set initial text values from default slider positions
    freqValueSpan.textContent = frequencySlider.value;
    pluckinessValueSpan.textContent = pluckinessSlider.value;
    modalPresetASelect.value = 'bell'; // Default Preset A
    modalPresetBSelect.value = 'levitationChamber';  // Default Preset B
    modalMorphSlider.value = 0.0;
    modalMorphValueSpan.textContent = parseFloat(modalMorphSlider.value).toFixed(2);
    chaosFeedbackSlider.value = 0.0;
    chaosFeedbackValueSpan.textContent = parseFloat(chaosFeedbackSlider.value).toFixed(2);
    delayTimeValueSpan.textContent = delayTimeSlider.value;
    delayFeedbackValueSpan.textContent = delayFeedbackSlider.value;
    filterCutoffValueSpan.textContent = filterCutoffSlider.value;
    masterGainValueSpan.textContent = masterGainSlider.value;
    mouseWheelControlTargetSelect.value = 'frequency'; // Default wheel target

    // Init NEW controls
    panLfoRateValueSpan.textContent = parseFloat(panLfoRateSlider.value).toFixed(2);
    panLfoDepthValueSpan.textContent = parseFloat(panLfoDepthSlider.value).toFixed(2);
    grainAmountValueSpan.textContent = parseFloat(grainAmountSlider.value).toFixed(3);

    // Set initial engine & update UI accordingly
    currentEngine = 'string';
    engineStringRadio.checked = true;
    updateUIForEngine(); // Sets damping label/range/value display + hides/shows controls

    // Set initial audio param values explicitly from sliders
    const now = audioContext.currentTime; // Use current time for scheduling initial values
    masterFilter.frequency.setValueAtTime(parseFloat(filterCutoffSlider.value), now);
    delay.delayTime.value = parseFloat(delayTimeSlider.value);
    delayFeedback.gain.setValueAtTime(parseFloat(delayFeedbackSlider.value), now);
    masterVolume.gain.setValueAtTime(parseFloat(masterGainSlider.value), now);
    // LFO params are set during their creation/initialization section above

    // Correctly assign engine-specific classes to parent control groups
    // Find the closest ancestor '.control-group' for robustness
    document.getElementById('pluckiness')?.closest('.control-group')?.classList.add('engine-specific', 'string-specific', 'hybrid-specific');
    modalSpecificDiv?.classList.add('engine-specific', 'modal-specific', 'hybrid-specific');
    // Damping and Frequency groups should remain universally visible (no engine-specific class)

    // Start the visualization loop!
    if (oscilloscopeCtx || spectrumCtx) {
        visualize();
     } else {
        console.warn("Canvas contexts not available, visualizations disabled.");
     }

}

// --- Run Initialization ---
initializeControlsAndAudio();
console.log("Advanced Synth Initialized. Features: Levitation Chamber, Pan LFO, Grain Delay. Play with keyboard or trigger button.");
console.log("Requires user interaction (click/key) to start audio.");