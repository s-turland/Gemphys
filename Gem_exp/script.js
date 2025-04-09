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

// Effect Control References
const panLfoRateSlider = document.getElementById('panLfoRate');
const panLfoDepthSlider = document.getElementById('panLfoDepth');
const grainAmountSlider = document.getElementById('grainAmount');

// Value display spans
const freqValueSpan = document.getElementById('freq-value');
const dampingValueSpan = document.getElementById('damping-value');
const pluckinessValueSpan = document.getElementById('pluckiness-value');
const modalMorphValueSpan = document.getElementById('modalMorph-value');
const chaosFeedbackValueSpan = document.getElementById('chaosFeedback-value');
const delayTimeValueSpan = document.getElementById('delayTime-value');
const delayFeedbackValueSpan = document.getElementById('delayFeedback-value');
const filterCutoffValueSpan = document.getElementById('filterCutoff-value');
const masterGainValueSpan = document.getElementById('masterGain-value');
const panLfoRateValueSpan = document.getElementById('panLfoRate-value');
const panLfoDepthValueSpan = document.getElementById('panLfoDepth-value');
const grainAmountValueSpan = document.getElementById('grainAmount-value');

// Visualizer References
const oscilloscopeCanvas = document.getElementById('oscilloscopeCanvas');
const spectrumCanvas = document.getElementById('spectrumCanvas');
let oscilloscopeCtx, spectrumCtx;

// MIDI References (Declared here, used within listeners/functions)
const midiDeviceSelector = document.getElementById('midiInDeviceSelector');
const midiStatus = document.getElementById('midiStatus');
let midiAccess = null; // Will hold the MIDIAccess object
let currentMIDIInput = null; // Will hold the currently selected MIDIInput device

// Mouse Wheel Control Target
const mouseWheelControlTargetSelect = document.getElementById('mouseWheelControlTarget');

// --- Polyphony and State Management ---
const MAX_VOICES = 10; // Allow slightly more for MIDI potentially
const activeVoices = []; // { key: midiNoteNumber | character, type, nodes, startTime, pannerNode, outputNode, gainNode } <- Added gainNode reference
const pressedKeys = new Set(); // Tracks pressed keys (both computer & MIDI note numbers)

// --- Keyboard Mapping & Tuning ---
const REFERENCE_MIDI_NOTE = 69; // A4 = MIDI note 69
const keyToNoteMap = {
     'a': -12, 's': -10, 'd': -9, 'f': -7, 'g': -5, 'h': -4, 'j': -2, 'k': 0, 'l': 2, ';': 3,
    'w': -11, 'e': -8,  't': -6, 'y': -3, 'u': -1, 'o': 1,  'p': 4
};

// Gets frequency relative to the A4 base frequency slider (MIDI note 69)
function getFrequencyFromMIDINote(midiNote) {
    const referenceFrequency = parseFloat(frequencySlider.value); // Use slider for A4 base
    // F = F_ref * 2^((m - m_ref)/12)
    return referenceFrequency * Math.pow(2, (midiNote - REFERENCE_MIDI_NOTE) / 12);
}

// Gets frequency for COMPUTER keyboard based on the slider freq for 'k' (MIDI 69 equivalent)
function getFrequencyFromKey(key) {
    const keyLower = key.toLowerCase();
    const semitoneDifference = keyToNoteMap[keyLower];
    if (semitoneDifference === undefined) return null;
    // Use the same MIDI note to freq calculation, relative to reference key 'k' (like A4=69)
    const referenceFrequency = parseFloat(frequencySlider.value);
    return referenceFrequency * Math.pow(2, semitoneDifference / 12);
}

// --- Modal Synthesis Data ---
const modalPresets = {
    bell: { modes: [ { ratio: 0.56, gain: 0.4, decay: 1.8 }, { ratio: 0.92, gain: 0.25, decay: 2.5 }, { ratio: 1.19, gain: 0.15, decay: 2.2 }, { ratio: 1.71, gain: 0.1, decay: 1.5 }, { ratio: 2.00, gain: 0.3, decay: 1.9 }, { ratio: 2.74, gain: 0.08, decay: 1.2 }, { ratio: 3.00, gain: 0.07, decay: 1.4 }, { ratio: 3.76, gain: 0.05, decay: 0.9 }, ], decayScale: 4.0 },
    bar: { modes: [ { ratio: 1.00, gain: 0.6, decay: 2.0 }, { ratio: 2.76, gain: 0.2, decay: 1.5 }, { ratio: 5.40, gain: 0.1, decay: 1.0 }, { ratio: 8.93, gain: 0.05, decay: 0.7 }, { ratio: 13.34, gain: 0.03, decay: 0.5 }, ], decayScale: 3.0 },
    glass: { modes: [ { ratio: 1.00, gain: 0.5, decay: 1.5 }, { ratio: 1.61, gain: 0.3, decay: 1.8 }, { ratio: 2.41, gain: 0.15, decay: 1.2 }, { ratio: 3.88, gain: 0.1, decay: 0.9 }, { ratio: 5.11, gain: 0.08, decay: 0.7 }, { ratio: 7.23, gain: 0.05, decay: 0.5 }, ], decayScale: 2.5 },
    tube: { modes: [ { ratio: 1.00, gain: 0.6, decay: 1.8 }, { ratio: 2.00, gain: 0.3, decay: 1.5 }, { ratio: 3.00, gain: 0.15, decay: 1.2 }, { ratio: 4.00, gain: 0.1, decay: 0.9 }, { ratio: 5.00, gain: 0.05, decay: 0.6 }, ], decayScale: 3.5 },
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

// --- Audio Routing Setup ---
const summingBus = audioContext.createGain();
const masterVolume = audioContext.createGain();
const masterFilter = audioContext.createBiquadFilter();
const delay = audioContext.createDelay(1.0); // Max delay time
const delayFeedback = audioContext.createGain();
const delayWetGain = audioContext.createGain(); // Controls overall delay level
const delayInputGain = audioContext.createGain(); // Input level to delay line

// LFOs and Controllers
const panLfo = audioContext.createOscillator();
const panLfoDepth = audioContext.createGain();
const delayLfo = audioContext.createOscillator();
const delayLfoDepth = audioContext.createGain(); // Controls grain modulation amount

// Analyser Node for Visualization
const analyserNode = audioContext.createAnalyser();
analyserNode.fftSize = 2048; // Standard FFT size
const timeDomainData = new Uint8Array(analyserNode.frequencyBinCount);
const frequencyData = new Uint8Array(analyserNode.frequencyBinCount);

// Connections
summingBus.connect(masterFilter);
masterFilter.connect(delayInputGain); // Feed post-filter signal to delay
masterFilter.connect(masterVolume); // Also connect post-filter signal to main output path

delayInputGain.connect(delay);
delay.connect(delayFeedback);
delayFeedback.connect(delay);
delay.connect(delayWetGain);
delayWetGain.connect(masterVolume); // Delay output rejoins main path AFTER filter, BEFORE master volume

masterVolume.connect(analyserNode); // Connect volume to analyser
analyserNode.connect(audioContext.destination); // Analyser connects to output

// LFO Connections (Modulators to Parameters)
panLfo.connect(panLfoDepth);
// panLfoDepth connects to each voice's pannerNode.pan parameter dynamically

delayLfo.connect(delayLfoDepth);
delayLfoDepth.connect(delay.delayTime); // Modulates the delay time

// Initial Effect Parameter Settings
masterFilter.type = 'lowpass';
masterFilter.Q.value = 1;
delayWetGain.gain.value = 0.7; // Overall delay mix
delayInputGain.gain.value = 0.6; // Signal fed into delay

// Start LFOs
panLfo.type = 'sine';
panLfo.frequency.value = parseFloat(panLfoRateSlider?.value ?? 0.5); // Use optional chaining and provide default
panLfoDepth.gain.value = parseFloat(panLfoDepthSlider?.value ?? 0);
panLfo.start();

delayLfo.type = 'triangle'; // Triangle or Sine often work well for smooth delay modulation
delayLfo.frequency.value = 6; // A typical rate for granular feel (adjust to taste)
delayLfoDepth.gain.value = parseFloat(grainAmountSlider?.value ?? 0);
delayLfo.start();


// --- Synthesis Functions (MODIFIED TO ACCEPT VELOCITY/RETURN GAIN NODE) ---

// Plucked String (Karplus-Strong) - Added velocity influence
function createPluckedStringVoice(frequency, damping, basePluckiness, velocity) { // Added velocity
    if (!audioContext || audioContext.state !== 'running') return null;
    const sampleRate = audioContext.sampleRate;
    const N = Math.round(sampleRate / frequency);
    const duration = 3.5;
    // Ensure bufferSize is valid, minimum N samples
    const bufferSize = Math.max(N, Math.round(sampleRate * duration));
    if (N <= 0 || bufferSize <=0) { console.warn("Invalid N or bufferSize for string voice"); return null; }
    const buffer = audioContext.createBuffer(1, bufferSize, sampleRate);
    const data = buffer.getChannelData(0);
    const delayLine = new Float32Array(N);

    // Velocity affects brightness (higher velocity = brighter pluck)
    const velocityFactor = (velocity / 127.0) * 0.6 + 0.7; // Map 0..127 -> 0.7..1.3 approx
    const effectivePluckiness = Math.max(0.1, Math.min(0.9, basePluckiness * velocityFactor));

    // Noise burst generation
    for (let i = 0; i < N; i++) delayLine[i] = Math.random() * 2 - 1;

    // Simple low-pass filter on noise based on effective pluckiness
    const pluckFilterFactor = 1.0 - effectivePluckiness; // Use velocity-modulated pluckiness
    let prevNoise = 0;
    for (let i = 0; i < N; i++) { let currentSample = delayLine[i]; delayLine[i] = currentSample * effectivePluckiness + prevNoise * pluckFilterFactor; prevNoise = delayLine[i]; }

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

    const sourceNode = audioContext.createBufferSource(); sourceNode.buffer = buffer;
    const gainNode = audioContext.createGain();
    // Use velocity to set initial gain (map 0..127 -> 0.0..0.9)
    const initialGain = (velocity / 127.0) * 0.9;
    gainNode.gain.setValueAtTime(initialGain, audioContext.currentTime); // Velocity controls level

    sourceNode.connect(gainNode);
    return { type: 'string', outputNode: gainNode, stopNodes: [sourceNode], gainNode: gainNode }; // Return gainNode
}

// Modal Synthesis - Added velocity influence
function createModalVoice(frequency, presetAData, presetBData, morphValue, globalDampingValue, chaosFeedbackAmount, velocity) { // Added velocity
     if (!audioContext || audioContext.state !== 'running') return null;
     if (!presetAData || !presetBData) { console.warn("Missing preset data"); return null; }

     const sampleRate = audioContext.sampleRate;
     const now = audioContext.currentTime;
     const voiceOutputGain = audioContext.createGain();
     // Set overall voice gain based on velocity
     const initialGain = (velocity / 127.0) * 0.9; // Map 0..127 -> 0..0.9
     voiceOutputGain.gain.setValueAtTime(initialGain, now);


     const oscillators = []; const gainNodes = []; const feedbackGainNodes = []; let longestDecay = 0;
     // Interpolate decay scale and calculate multiplier based on global damping
     const interpDecayScale = presetAData.decayScale + (presetBData.decayScale - presetAData.decayScale) * morphValue;
     // Inverse relationship for modal: Lower damping slider value = longer decay
     const decayMultiplier = Math.max(0.1, interpDecayScale * (1.05 - globalDampingValue));
     const numModes = Math.max(presetAData.modes.length, presetBData.modes.length);

     for (let i = 0; i < numModes; i++) {
         // Interpolate mode parameters (ratio, gain, decay)
        const modeA = presetAData.modes[i];
        const modeB = presetBData.modes[i];
        let interpRatio, interpGain, interpDecay, effectiveGain;
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
        } else { continue; } // Skip if mode doesn't exist in either

        const modeFreq = frequency * interpRatio;
        const decayTime = interpDecay * decayMultiplier; // Apply calculated decay

        // Skip invalid or inaudible modes
        if (modeFreq > sampleRate / 2 || modeFreq <= 0 || effectiveGain < 0.001 || decayTime <= 0.001) continue;

        const osc = audioContext.createOscillator(); const gain = audioContext.createGain(); // This mode's primary gain envelope
        osc.frequency.setValueAtTime(modeFreq, now); osc.type = 'sine';

        // Gain envelope: quick attack, exponential decay
        // Set initial mode gain - keep effectiveGain as the main factor here, overall level set by voiceOutputGain
        gain.gain.setValueAtTime(0, now);
        gain.gain.linearRampToValueAtTime(effectiveGain, now + 0.005); // Short attack
        gain.gain.setTargetAtTime(0.0001, now + 0.01, decayTime / 4); // Exponential decay target

        osc.connect(gain); gain.connect(voiceOutputGain);
        osc.start(now); const stopTime = now + decayTime * 1.5 + 0.2; // Schedule stop well after decay finishes
        osc.stop(stopTime);
        if (decayTime > longestDecay) longestDecay = decayTime;
        oscillators.push(osc); gainNodes.push(gain); // Store for chaos feedback
     }

    // Chaotic Feedback Implementation
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
    }
    if (oscillators.length === 0) { console.warn("No valid modes created for modal voice"); return null;}

     return { type: 'modal', outputNode: voiceOutputGain, stopNodes: oscillators, feedbackNodes: feedbackGainNodes, gainNode: voiceOutputGain }; // Return main gain node
}

// Hybrid Resonance - Added velocity influence
function createHybridResonanceVoice(frequency, presetAData, presetBData, morphValue, baseExcitationPluckiness, modalDampingValue, chaosFeedbackAmount, velocity) { // Added velocity
    if (!audioContext || audioContext.state !== 'running') return null;
    if (!presetAData || !presetBData) { console.warn("Missing preset data for hybrid voice"); return null; }

    const sampleRate = audioContext.sampleRate; const now = audioContext.currentTime;
    const voiceOutputGain = audioContext.createGain();
    // Velocity sets overall gain
    const initialGain = (velocity / 127.0) * 0.9; // Map 0..127 -> 0..0.9
    voiceOutputGain.gain.setValueAtTime(initialGain, now);

    const oscillators = []; const primaryGainNodes = []; const feedbackGainNodes = []; let longestDecay = 0;

    // --- Virtual String Excitation Simulation ---
    // Velocity influences pluckiness
    const velocityFactor = (velocity / 127.0) * 0.6 + 0.7; // Map 0..127 -> 0.7..1.3
    const excitationPluckiness = Math.max(0.1, Math.min(0.9, baseExcitationPluckiness * velocityFactor));

    const virtualStringFundamental = frequency;
    const VIRTUAL_STRING_BASE_DECAY = 1.5; // Base decay of fundamental (seconds)
    const VIRTUAL_STRING_DAMPING_PER_OCTAVE = 0.5; // Higher harmonics decay faster
    const MAX_VIRTUAL_HARMONICS = 20;
    const virtualHarmonicAmplitudes = []; // Calculated amplitude for each harmonic
    const baseAmplitudes = []; let amplitudeSum = 0;

    // Initial amplitudes (1/n)
    for(let i = 1; i <= MAX_VIRTUAL_HARMONICS; i++) { baseAmplitudes[i] = 1 / i; }

    // Simulate pluckiness filtering on harmonic amplitudes (simple LPF using velocity-modulated pluckiness)
    const filterFactor = 1.0 - excitationPluckiness; // USE MODULATED VALUE
    let previousFilteredAmp = 0;
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
        // Interpolate Mode Parameters (Same logic as in createModalVoice)
        const modeA = presetAData.modes[i]; const modeB = presetBData.modes[i];
        let interpRatio, interpGain, interpDecay, effectiveBaseGain;
         if (modeA && modeB) {
              interpRatio = modeA.ratio + (modeB.ratio - modeA.ratio) * morphValue; interpGain = modeA.gain + (modeB.gain - modeA.gain) * morphValue; interpDecay = modeA.decay + (modeB.decay - modeA.decay) * morphValue; effectiveBaseGain = interpGain;
          } else if (modeA) { interpRatio = modeA.ratio; interpGain = modeA.gain; interpDecay = modeA.decay; effectiveBaseGain = interpGain * (1.0 - morphValue); }
          else if (modeB) { interpRatio = modeB.ratio; interpGain = modeB.gain; interpDecay = modeB.decay; effectiveBaseGain = interpGain * morphValue; }
          else { continue; }
          // --- End Interpolation ---

        const osc = audioContext.createOscillator();
        const primaryModeGain = audioContext.createGain(); // Controls the modal decay envelope
        const activationGain = audioContext.createGain(); // Controls amplitude based on string excitation

        const modeFreq = frequency * interpRatio;
        const baseModeDecayTime = interpDecay * decayMultiplier; // Base decay for this mode

        if (modeFreq > sampleRate / 2 || modeFreq <= 0 || effectiveBaseGain < 0.001 || baseModeDecayTime <= 0.001) continue;

        osc.frequency.setValueAtTime(modeFreq, now);
        osc.type = 'sine';

        // --- Calculate Virtual String Influence (uses velocity-modded amplitudes) ---
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
        // Initial gain is set by main voiceOutputGain based on velocity

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

    // Chaotic Feedback Implementation (Connects primaryGain nodes)
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
    }


    if (oscillators.length === 0) { console.warn(`No valid modes/harmonics for hybrid. Freq: ${frequency.toFixed(1)}, Morph:${morphValue.toFixed(2)}`); return null; }

    return { type: 'hybrid', outputNode: voiceOutputGain, stopNodes: oscillators, feedbackNodes: feedbackGainNodes, gainNode: voiceOutputGain }; // Return main gain node
}


// --- Note Triggering Functions (Adapted for both MIDI and Keyboard) ---

function playNote(keyIdentifier, frequency, velocity = 90) { // Added velocity, default reasonable value
     if (!isAudioContextResumed || audioContext.state !== 'running') {
        console.log("Audio Context not running, trying to resume for", keyIdentifier);
        if (audioContext.state === 'suspended') {
            resumeAudioContext().then(() => {
                if (audioContext.state === 'running') playNote(keyIdentifier, frequency, velocity);
            });
        }
        return;
     }
     if (!frequency || frequency <= 0) {
        console.warn("Attempted to play note with invalid frequency:", frequency, "for key:", keyIdentifier);
        return;
     }

     // --- Voice Management (Stealing & Retriggering) --- Check if keyIdentifier is already playing
    if (pressedKeys.has(keyIdentifier)) {
        // console.log("Retriggering note:", keyIdentifier);
        stopNote(keyIdentifier, 0.005); // Quick stop before restart
    }

    if (activeVoices.length >= MAX_VOICES) {
        activeVoices.sort((a, b) => a.startTime - b.startTime); // Find oldest voice
        const stolenVoiceData = activeVoices.shift(); // Remove it from list immediately
        console.log(`Max voices reached. Stealing oldest voice (Key: ${stolenVoiceData.key}, Type: ${stolenVoiceData.type}).`);
        stopNote(stolenVoiceData.key, 0.015); // Trigger fade out for the stolen voice's key ID
        // Note: We remove the voice from activeVoices *here* to make room. stopNote handles the actual audio stop and pressedKeys.
    }


    // --- Create Voice Based on Engine ---
    let voiceData = null;
    const now = audioContext.currentTime;
    const basePluck = parseFloat(pluckinessSlider?.value ?? 0.5); // Use optional chaining and default
    const dampingValue = parseFloat(dampingSlider?.value ?? (currentEngine === 'string' ? 0.996 : 0.5));

    if (currentEngine === 'string') {
        voiceData = createPluckedStringVoice(frequency, dampingValue, basePluck, velocity);
    } else { // Modal or Hybrid
        const presetAName = modalPresetASelect?.value ?? 'bell';
        const presetBName = modalPresetBSelect?.value ?? 'levitationChamber';
        const presetAData = modalPresets[presetAName];
        const presetBData = modalPresets[presetBName];
        const morphValue = parseFloat(modalMorphSlider?.value ?? 0);
        const chaosValue = parseFloat(chaosFeedbackSlider?.value ?? 0);

        if (!presetAData || !presetBData) {
             console.error("Could not find preset data for morphing!");
             return; // Stop if presets are missing
        }

        if (currentEngine === 'modal') {
            voiceData = createModalVoice(frequency, presetAData, presetBData, morphValue, dampingValue, chaosValue, velocity);
        } else { // Hybrid engine
            voiceData = createHybridResonanceVoice(frequency, presetAData, presetBData, morphValue, basePluck, dampingValue, chaosValue, velocity);
        }
    }


    // --- Connect and Store Voice ---
    if (voiceData) {
        const pannerNode = audioContext.createStereoPanner();
        panLfoDepth.connect(pannerNode.pan); // Connect the GLOBAL Pan LFO Depth controller to THIS voice's panner parameter
        voiceData.outputNode.connect(pannerNode);
        pannerNode.connect(summingBus); // Voice Output -> Panner -> Summing Bus

        if (voiceData.type === 'string' && voiceData.stopNodes && voiceData.stopNodes[0]) {
             voiceData.stopNodes[0].start(now); // Start buffer source immediately
             // Cleanup for buffer sources needs onended
             voiceData.stopNodes[0].onended = () => cleanupStoppedVoice(voiceData.stopNodes[0]);
        }
        // Modal & Hybrid oscillators rely on scheduled stop() or explicit note off/stealing for cleanup

        const newVoice = {
             key: keyIdentifier, // Use MIDI note number or character
             type: voiceData.type,
             nodes: voiceData,
             outputNode: voiceData.outputNode, // The final gain node for the voice
             pannerNode: pannerNode,           // Panner specific to this voice
             gainNode: voiceData.gainNode,     // Reference to the gain node velocity affected
             startTime: now
         };
         activeVoices.push(newVoice);
         pressedKeys.add(keyIdentifier); // Track pressed state

     } else {
          console.log("Failed to create voice data for key:", keyIdentifier);
     }
}

// Function to explicitly stop a note (e.g., on MIDI Note Off or voice stealing)
function stopNote(keyIdentifier, releaseTime = 0.05) { // Add releaseTime parameter
    // Find the *last* started voice for this key, as that's most likely the one actively sounding.
    let voiceIndex = -1;
    let voiceToStop = null;
    for (let i = activeVoices.length - 1; i >= 0; i--) {
        if (activeVoices[i].key === keyIdentifier) {
            voiceIndex = i;
            voiceToStop = activeVoices[i];
            break;
        }
    }

    // Remove from pressed keys immediately, regardless of finding an active voice
    // This allows retriggering even if cleanup is slightly delayed or fails.
    pressedKeys.delete(keyIdentifier);

    if (voiceToStop) {
        // console.log(`Stopping note: ${keyIdentifier}, Voice Type: ${voiceToStop.type}, StartTime: ${voiceToStop.startTime}`);

        // --- Smooth Fade Out using Gain Node ---
        if (voiceToStop.gainNode && voiceToStop.gainNode.gain) {
            const now = audioContext.currentTime;
             try {
                 // Check if gain param exists before cancelling/scheduling
                 if(voiceToStop.gainNode.gain) {
                    voiceToStop.gainNode.gain.cancelScheduledValues(now);
                    voiceToStop.gainNode.gain.setTargetAtTime(0.0, now, releaseTime / 3); // Faster exponential fade
                 } else {
                    console.warn("Gain node parameter missing for stopNote:", keyIdentifier);
                 }
             } catch (e) {
                 console.warn("Error scheduling gain for stopNote:", keyIdentifier, e);
                 // Node might have been disconnected already by stealing mechanism etc.
             }
        } else {
            console.warn("Could not find gain node for fade out on voice:", keyIdentifier);
        }

        // Schedule cleanup slightly after fade to ensure audio stops before disconnection
        const cleanupDelay = releaseTime + 0.15; // Increased buffer time slightly
        setTimeout(() => {
            // Find the *exact same* voice instance again using startTime for safety before removing & disconnecting
            const currentIndex = activeVoices.findIndex(v => v.startTime === voiceToStop.startTime); // Use startTime as unique ID
             if (currentIndex > -1) {
                const stoppedVoice = activeVoices.splice(currentIndex, 1)[0]; // Remove from active list
                 // console.log(`Cleaning up stopped note: ${stoppedVoice.key}`);
                 try { panLfoDepth.disconnect(stoppedVoice.pannerNode.pan); } catch(e) {}
                 try { stoppedVoice.pannerNode.disconnect(); } catch(e) {}
                 try { stoppedVoice.outputNode.disconnect(); } catch(e) {} // Disconnect main output too

                 // Stop oscillators/sources if they haven't finished their decay/scheduled stop
                 stoppedVoice.nodes.stopNodes.forEach(node => {
                    if (node && typeof node.stop === 'function') {
                        try {
                            // Don't stop immediately, let fade complete. Stop slightly in future.
                            node.stop(audioContext.currentTime + 0.05);
                            node.disconnect();
                        } catch(e){
                            // Ignore errors if already stopped/disconnected
                        }
                    }
                });
                // Disconnect feedback nodes too if they exist
                 if(stoppedVoice.nodes.feedbackNodes) {
                    stoppedVoice.nodes.feedbackNodes.forEach(node => { try { node.disconnect(); } catch(e) {} });
                 }

             } else {
                // console.log("Voice to cleanup already removed or mismatch:", keyIdentifier, voiceToStop.startTime);
             }
        }, cleanupDelay * 1000); // setTimeout expects milliseconds

    } else {
         // console.log("Stop request for inactive note (or voice already stolen/cleaned up):", keyIdentifier);
    }
}

// Cleanup function specifically for BufferSource onended (String engine)
function cleanupStoppedVoice(stoppedNode) {
    // Find the voice associated with the stopped buffer source
    const voiceIndex = activeVoices.findIndex(v => v.nodes && v.nodes.stopNodes && v.nodes.stopNodes.includes(stoppedNode));
    if (voiceIndex > -1) {
        const removedVoice = activeVoices.splice(voiceIndex, 1)[0]; // Remove from active list
        // console.log("Cleaning up ended buffer source voice:", removedVoice.key);
        // Disconnect nodes
        try { panLfoDepth.disconnect(removedVoice.pannerNode.pan); } catch(e) {}
        try { removedVoice.pannerNode.disconnect(); } catch(e) {}
        try { removedVoice.outputNode.disconnect(); } catch(e) {}
        // Ensure it's removed from pressedKeys if it ended naturally without a note-off
        pressedKeys.delete(removedVoice.key);
    }
}


// --- MIDI Handling ---

function onMIDISuccess(mAccess) {
    console.log("MIDI ready!");
    midiAccess = mAccess;
    listMIDIDevices();
    midiAccess.onstatechange = listMIDIDevices; // Update list if devices change/added/removed
}

function onMIDIFailure(msg) {
    console.error(`Failed to get MIDI access - ${msg}`);
    // Use the globally declared midiStatus element
    if (midiStatus) {
        midiStatus.textContent = "MIDI Error: " + msg;
        midiStatus.style.color = 'var(--accent-color-warning)';
    }
    if (midiDeviceSelector) {
        midiDeviceSelector.disabled = true;
    }
}

function listMIDIDevices() {
    if (!midiAccess || !midiDeviceSelector || !midiStatus) {
        console.warn("MIDI Access or UI elements not available for listing devices.");
        return;
    }

    // Clear existing options (except the default "-- Select --")
    while (midiDeviceSelector.options.length > 1) {
        midiDeviceSelector.remove(1);
    }

    let count = 0;
    if (midiAccess.inputs.size > 0) {
        midiAccess.inputs.forEach((input) => {
            const option = document.createElement("option");
            option.value = input.id;
            option.text = input.name;
            midiDeviceSelector.appendChild(option);
            count++;
        });
         midiStatus.textContent = `${count} MIDI device(s) found. Select one.`;
         midiStatus.style.color = 'var(--label-color)';
         midiDeviceSelector.disabled = false;

         // Restore previously selected device if it still exists
         if (currentMIDIInput && midiAccess.inputs.has(currentMIDIInput.id)) {
              midiDeviceSelector.value = currentMIDIInput.id;
              // Status will be updated by connectMIDIInput if called again, otherwise show connected status
              if (currentMIDIInput.connection === "open") { // Check connection state
                    midiStatus.textContent = `Connected: ${currentMIDIInput.name}`;
                    midiStatus.style.color = 'var(--accent-color-positive)';
              } else {
                    // It exists but isn't connected? Maybe prompt selection again.
                    midiStatus.textContent = "Select a MIDI device.";
                    midiDeviceSelector.value = "";
              }
         } else {
             midiDeviceSelector.value = ""; // Reset selection
             disconnectCurrentMIDI(); // Disconnect if previous device gone or wasn't selected
             midiStatus.textContent = "Select a MIDI device."; // Prompt user
         }

    } else {
        midiStatus.textContent = "No MIDI input devices found.";
        midiStatus.style.color = 'var(--accent-color-warning)';
        midiDeviceSelector.disabled = true;
        disconnectCurrentMIDI();
    }
}

function connectMIDIInput(deviceId) {
    if (!midiAccess || !midiStatus) return; // Need access and status element

    // Handle case where user selects "-- Select --"
    if (!deviceId) {
        disconnectCurrentMIDI();
        midiStatus.textContent = "Select a MIDI device.";
        midiStatus.style.color = 'var(--label-color)';
        if (midiDeviceSelector) midiDeviceSelector.value = "";
        return;
    }

    disconnectCurrentMIDI(); // Disconnect previous listener first

    const selectedInput = midiAccess.inputs.get(deviceId);
    if (selectedInput) {
        console.log(`Connecting to MIDI device: ${selectedInput.name} (ID: ${selectedInput.id})`);

        // Add event listener and open the port
        selectedInput.onmidimessage = handleMIDIMessage;
        selectedInput.open()
            .then(() => {
                currentMIDIInput = selectedInput; // Store reference *after* successful open
                midiStatus.textContent = `Connected: ${selectedInput.name}`;
                midiStatus.style.color = 'var(--accent-color-positive)';
                console.log(`Successfully opened and connected to ${selectedInput.name}`);
            })
            .catch((err) => {
                 console.error(`Failed to open MIDI port ${selectedInput.name}:`, err);
                 midiStatus.textContent = `Error connecting to ${selectedInput.name}.`;
                 midiStatus.style.color = 'var(--accent-color-warning)';
                 selectedInput.onmidimessage = null; // Remove listener if open failed
                 currentMIDIInput = null;
                 if (midiDeviceSelector) midiDeviceSelector.value = ""; // Reset dropdown
            });

    } else {
        console.warn("Could not find selected MIDI device with ID:", deviceId);
         midiStatus.textContent = "Could not connect to device.";
         midiStatus.style.color = 'var(--accent-color-warning)';
         currentMIDIInput = null; // Clear reference
         if (midiDeviceSelector) midiDeviceSelector.value = ""; // Reset dropdown
    }
}

function disconnectCurrentMIDI() {
     if (currentMIDIInput) {
        console.log(`Disconnecting from MIDI device: ${currentMIDIInput.name}`);
        currentMIDIInput.onmidimessage = null; // Remove listener
        currentMIDIInput.close() // Close the port (returns a promise, but we don't necessarily wait)
             .catch(e => console.warn("Error closing MIDI port:", e));
        currentMIDIInput = null;
        // Status will be updated by listMIDIDevices or connectMIDIInput
     }
}

function handleMIDIMessage(event) {
    // Ensure audio context is running before processing MIDI
    if (!isAudioContextResumed || audioContext.state !== 'running') {
        // Optionally try to resume, but might flood console if many messages arrive while suspended
        // resumeAudioContext(); // Consider implications before enabling this here
        return; // Don't process MIDI if context isn't ready
    }

    const command = event.data[0] >> 4; // Command type (ignore channel for now)
    const channel = event.data[0] & 0xf; // Channel (0-15) - unused currently
    const note = event.data[1];        // MIDI note number (0-127)
    const velocity = event.data.length > 2 ? event.data[2] : 0; // Velocity (0-127)

    // Note On (Command 9) - Check for velocity > 0
    if (command === 9 && velocity > 0) {
        const frequency = getFrequencyFromMIDINote(note);
        if (!frequency) return; // Avoid playing if frequency is invalid
        // console.log(`MIDI Note On: Note=${note}, Vel=${velocity}, Freq=${frequency.toFixed(2)}`);
        playNote(note, frequency, velocity); // Use MIDI note number as keyIdentifier
    }
    // Note Off (Command 8) OR Note On with Velocity 0 (common practice)
    else if (command === 8 || (command === 9 && velocity === 0)) {
         // console.log(`MIDI Note Off: Note=${note}`);
         stopNote(note); // Use MIDI note number as keyIdentifier to stop
    }
    // --- Control Change (CC) Handling ---
    else if (command === 11) { // Control Change (0xB0-0xBF)
        const ccNumber = event.data[1];
        const ccValue = event.data[2]; // 0-127
        // console.log(`MIDI CC: Ch=${channel} CC=${ccNumber} Val=${ccValue}`);
        handleMIDIControlChange(ccNumber, ccValue); // Delegate to separate function
    }
}

function handleMIDIControlChange(ccNumber, ccValue) {
     // Map 0-127 value to 0.0-1.0 range for easier scaling
     const valueNormalized = ccValue / 127.0;

     // Find the corresponding slider/control based on CC number
     let targetSlider = null;
     let targetParamSetter = null; // Optional direct audio param function

     switch(ccNumber) {
        case 1: // Mod Wheel - Map to Morph
             targetSlider = modalMorphSlider;
             break;
        case 74: // Filter Cutoff (Common assignment)
             targetSlider = filterCutoffSlider;
             break;
        case 71: // Filter Resonance/Q (Common assignment)
             // Direct mapping example (assuming Q range 0.1 to 15)
             const qValue = 0.1 + valueNormalized * 14.9;
             targetParamSetter = () => masterFilter.Q.setTargetAtTime(qValue, audioContext.currentTime, 0.01);
             // console.log(`MIDI CC ${ccNumber} -> Filter Q: ${qValue.toFixed(2)}`);
             break;
        case 10: // Pan (Common assignment) - Map to Pan LFO Depth maybe?
             targetSlider = panLfoDepthSlider;
             break;
        case 7: // Volume (Common assignment) - Map to Master Gain
             targetSlider = masterGainSlider;
             break;
        case 11: // Expression (Common assignment) - Also map to Master Gain? Or Filter?
             // Example: Map expression to master gain as well
             // targetSlider = masterGainSlider;
             break;
        // Add more CC mappings here...
        // case 80: // Example: Map to Chaos Feedback
        //     targetSlider = chaosFeedbackSlider;
        //     break;
        // case 81: // Example: Map to Grain Amount
        //     targetSlider = grainAmountSlider;
        //     break;
     }

     // If a target slider was found, update its value
     if (targetSlider) {
         const min = parseFloat(targetSlider.min);
         const max = parseFloat(targetSlider.max);
         const newValue = min + valueNormalized * (max - min);

         // Check if slider exists and is visible before updating
         const sliderGroup = targetSlider.closest('.control-group');
         const style = window.getComputedStyle(sliderGroup || targetSlider);
         if (style && style.display !== 'none' && style.visibility !== 'hidden') {
            targetSlider.value = newValue;
            targetSlider.dispatchEvent(new Event('input', { bubbles: true })); // Trigger UI/Audio update
            // console.log(`MIDI CC ${ccNumber} -> ${targetSlider.id}: ${newValue.toFixed(2)}`);
         } else {
            // console.log(`MIDI CC ${ccNumber}: Target slider ${targetSlider.id} is hidden, ignoring.`);
         }
     }
     // If a direct parameter setter exists, call it
     else if (targetParamSetter) {
         targetParamSetter();
     }
}


// --- Visualization Functions ---
function drawOscilloscope() {
    if (!oscilloscopeCtx || !analyserNode || !oscilloscopeCanvas) return;
    const width = oscilloscopeCanvas.clientWidth; // Use clientWidth for responsive sizing
    const height = oscilloscopeCanvas.clientHeight;

    // Check if canvas has valid dimensions
    if (width <= 0 || height <= 0) return;

    // Ensure internal buffer matches client size * dpr (done in setupCanvas)
    // oscilloscopeCtx.setTransform(window.devicePixelRatio || 1, 0, 0, window.devicePixelRatio || 1, 0, 0); // Set in setupCanvas

    oscilloscopeCtx.clearRect(0, 0, width, height); // Clear using client dimensions scaled by ctx

    analyserNode.getByteTimeDomainData(timeDomainData);

    oscilloscopeCtx.lineWidth = 1.5;
    oscilloscopeCtx.strokeStyle = 'rgb(200, 255, 200)'; // Light green color
    oscilloscopeCtx.beginPath();

    const bufferLength = analyserNode.frequencyBinCount; // This is actually half fftSize
    const sliceWidth = width / bufferLength;
    let x = 0;

    for (let i = 0; i < bufferLength; i++) {
        const v = timeDomainData[i] / 128.0; // Normalize to 0-2 range
        const y = v * height / 2;

        if (i === 0) {
            oscilloscopeCtx.moveTo(x, y);
        } else {
            oscilloscopeCtx.lineTo(x, y);
        }
        x += sliceWidth;
    }

    oscilloscopeCtx.lineTo(width, height / 2); // End line at center height
    oscilloscopeCtx.stroke();
}

function drawSpectrum() {
    if (!spectrumCtx || !analyserNode || !spectrumCanvas) return;
    const width = spectrumCanvas.clientWidth;
    const height = spectrumCanvas.clientHeight;

    if (width <= 0 || height <= 0) return;

    spectrumCtx.clearRect(0, 0, width, height);

    analyserNode.getByteFrequencyData(frequencyData);

    const bufferLength = analyserNode.frequencyBinCount; // number of data points
    const barWidth = (width / bufferLength) * 1.5; // Slightly wider bars
    let barHeight;
    let x = 0;

    // Use a gradient for visual appeal
    const gradient = spectrumCtx.createLinearGradient(0, height, 0, 0);
    gradient.addColorStop(0, "rgb(50, 150, 50)");   // Darker Green at bottom
    gradient.addColorStop(0.6, "rgb(100, 200, 100)"); // Mid Green
    gradient.addColorStop(1, "rgb(180, 255, 180)"); // Lighter Green at top
    spectrumCtx.fillStyle = gradient; // Use gradient

    for (let i = 0; i < bufferLength; i++) {
        barHeight = frequencyData[i] * (height / 255.0); // Scale bar height

        spectrumCtx.fillRect(x, height - barHeight, barWidth, barHeight);

        x += barWidth + 1; // Add spacing between bars
    }
}

let visualizationFrameId = null;
function visualize() {
    if (oscilloscopeCtx) drawOscilloscope();
    if (spectrumCtx) drawSpectrum();
    visualizationFrameId = requestAnimationFrame(visualize); // Loop the visualization
}
function stopVisualize() {
    if(visualizationFrameId) {
        cancelAnimationFrame(visualizationFrameId);
        visualizationFrameId = null;
    }
}


// --- Function to Resume Audio Context ---
async function resumeAudioContext() {
    if (audioContext.state === 'suspended') {
        console.log('Attempting to resume AudioContext...');
        try {
            await audioContext.resume();
            if (audioContext.state === 'running') {
                 isAudioContextResumed = true;
                 console.log('AudioContext is running.');
                 // Restart LFOs if they were stopped or not started properly
                 // Note: LFOs are usually started once and run continuously.
            } else {
                isAudioContextResumed = false;
                console.warn('AudioContext could not be resumed. State:', audioContext.state);
            }
        } catch (e) {
            console.error("Error resuming AudioContext:", e);
            isAudioContextResumed = false;
        }
    } else if (audioContext.state === 'running') {
        isAudioContextResumed = true; // Already running
    } else {
        isAudioContextResumed = false; // Other states like 'closed'
    }
}

// --- Update UI based on Engine --- Handles dynamic labels/ranges/visibility ---
function updateUIForEngine() {
    const body = document.body;

    // --- Step 1: Clear existing engine classes ---
    body.classList.remove('engine-string', 'engine-modal', 'engine-hybrid');

    // --- Step 2: Get references (ensure elements exist) ---
    const modalHybridWheelOptions = document.querySelectorAll('#mouseWheelControlTarget .modal-hybrid-wheel');

    // Check essential elements exist before proceeding
    if (!dampingLabel || !dampingHelp || !pluckinessLabel || !pluckinessHelp || !dampingSlider || !pluckinessSlider || !dampingValueSpan || !pluckinessValueSpan) {
        console.error("Essential UI elements for engine update are missing!");
        return;
    }

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
        let currentDampVal = parseFloat(dampingSlider.value);
        if (isNaN(currentDampVal) || currentDampVal < 0.95 || currentDampVal > 0.999) {
            dampingSlider.value = 0.996; // Default string damping
        }

        // Update Mouse Wheel Options Visibility
        modalHybridWheelOptions.forEach(opt => opt.style.display = 'none');
        // Switch wheel target if current selection is now hidden
        if (mouseWheelControlTargetSelect && ['modalMorph', 'chaosFeedback'].includes(mouseWheelControlTargetSelect.value)) {
            mouseWheelControlTargetSelect.value = 'frequency'; // Default safe choice
        }

    } else { // Modal or Hybrid (Share many UI properties)
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
        let currentDampVal = parseFloat(dampingSlider.value);
        if (isNaN(currentDampVal) || currentDampVal < 0.1 || currentDampVal > 1.0) {
            dampingSlider.value = 0.5; // Default mid-value for modal/hybrid
        }

        // Update Mouse Wheel Options Visibility
         modalHybridWheelOptions.forEach(opt => opt.style.display = 'block'); // Show morph/chaos
    }

    // --- Step 4: Update value displays AFTER potentially changing slider values ---
    const dampVal = parseFloat(dampingSlider.value);
    dampingValueSpan.textContent = (currentEngine === 'string' ? dampVal.toFixed(3) : dampVal.toFixed(2));

    // Update pluckiness display only if its control is relevant
    if (currentEngine === 'string' || currentEngine === 'hybrid') {
        pluckinessValueSpan.textContent = parseFloat(pluckinessSlider.value).toFixed(2);
    } else {
        pluckinessValueSpan.textContent = "---"; // Indicate N/A when hidden
    }

    // --- Step 5: Rely on CSS rules linked to body classes ---
    // This part is handled by the CSS based on body.engine-* classes.
}


// --- Event Listeners Setup Function ---
// Encapsulate listener setup to be called after DOM loaded
function setupEventListeners() {
    // Engine Switch
    if(engineStringRadio) engineStringRadio.addEventListener('change', () => { if (engineStringRadio.checked) { currentEngine = 'string'; updateUIForEngine(); } });
    if(engineModalRadio) engineModalRadio.addEventListener('change', () => { if (engineModalRadio.checked) { currentEngine = 'modal'; updateUIForEngine(); } });
    if(engineHybridRadio) engineHybridRadio.addEventListener('change', () => { if (engineHybridRadio.checked) { currentEngine = 'hybrid'; updateUIForEngine(); } });

    // Pluck Button
    if(pluckButton) pluckButton.addEventListener('click', async () => {
        await resumeAudioContext(); // Ensure context is running
        if (isAudioContextResumed) {
            const referenceFrequency = parseFloat(frequencySlider?.value ?? 440);
            // Use the MIDI reference note for consistency, trigger with default velocity
            playNote(REFERENCE_MIDI_NOTE, referenceFrequency, 90);
        }
    });

    // Slider Updates (Connect UI to Audio Params) - Check elements exist
    if(frequencySlider && freqValueSpan) frequencySlider.oninput = () => freqValueSpan.textContent = frequencySlider.value;
    if(dampingSlider && dampingValueSpan) dampingSlider.oninput = () => {
        const dampVal = parseFloat(dampingSlider.value);
        dampingValueSpan.textContent = (currentEngine === 'string' ? dampVal.toFixed(3) : dampVal.toFixed(2));
    };
    if(pluckinessSlider && pluckinessValueSpan) pluckinessSlider.oninput = () => pluckinessValueSpan.textContent = parseFloat(pluckinessSlider.value).toFixed(2);
    if(masterGainSlider && masterGainValueSpan) masterGainSlider.oninput = () => { masterGainValueSpan.textContent = masterGainSlider.value; masterVolume.gain.setTargetAtTime(parseFloat(masterGainSlider.value), audioContext.currentTime, 0.01); };
    if(delayTimeSlider && delayTimeValueSpan) delayTimeSlider.oninput = () => { delayTimeValueSpan.textContent = delayTimeSlider.value; delay.delayTime.setTargetAtTime(parseFloat(delayTimeSlider.value), audioContext.currentTime, 0.01); };
    if(delayFeedbackSlider && delayFeedbackValueSpan) delayFeedbackSlider.oninput = () => { delayFeedbackValueSpan.textContent = delayFeedbackSlider.value; delayFeedback.gain.setTargetAtTime(parseFloat(delayFeedbackSlider.value), audioContext.currentTime, 0.01); };
    if(filterCutoffSlider && filterCutoffValueSpan) filterCutoffSlider.oninput = () => { filterCutoffValueSpan.textContent = filterCutoffSlider.value; masterFilter.frequency.setTargetAtTime(parseFloat(filterCutoffSlider.value), audioContext.currentTime, 0.01); };
    if(modalMorphSlider && modalMorphValueSpan) modalMorphSlider.oninput = () => modalMorphValueSpan.textContent = parseFloat(modalMorphSlider.value).toFixed(2);
    if(chaosFeedbackSlider && chaosFeedbackValueSpan) chaosFeedbackSlider.oninput = () => chaosFeedbackValueSpan.textContent = parseFloat(chaosFeedbackSlider.value).toFixed(2);
    if(panLfoRateSlider && panLfoRateValueSpan) panLfoRateSlider.oninput = () => { panLfoRateValueSpan.textContent = parseFloat(panLfoRateSlider.value).toFixed(2); panLfo.frequency.setTargetAtTime(parseFloat(panLfoRateSlider.value), audioContext.currentTime, 0.01); };
    if(panLfoDepthSlider && panLfoDepthValueSpan) panLfoDepthSlider.oninput = () => { panLfoDepthValueSpan.textContent = parseFloat(panLfoDepthSlider.value).toFixed(2); panLfoDepth.gain.setTargetAtTime(parseFloat(panLfoDepthSlider.value), audioContext.currentTime, 0.01); };
    if(grainAmountSlider && grainAmountValueSpan) grainAmountSlider.oninput = () => { grainAmountValueSpan.textContent = parseFloat(grainAmountSlider.value).toFixed(3); delayLfoDepth.gain.setTargetAtTime(parseFloat(grainAmountSlider.value), audioContext.currentTime, 0.01); };

    // Mouse Wheel Control
    if(controlsDiv && mouseWheelControlTargetSelect) {
         controlsDiv.addEventListener('wheel', handleWheelControl, { passive: false }); // passive: false needed for preventDefault
    }

    // Add a general click/key listener to resume context if needed (safer than adding on body initially)
    document.body.addEventListener('click', resumeAudioContext, { once: true });
    document.body.addEventListener('keydown', resumeAudioContext, { once: true });

    // Window resize listener for canvases
    window.addEventListener('resize', () => {
        setupCanvas(oscilloscopeCanvas, oscilloscopeCtx);
        setupCanvas(spectrumCanvas, spectrumCtx);
    });
}

// --- Mouse Wheel Control ---
function handleWheelControl(event) {
    if (!mouseWheelControlTargetSelect) return; // Target selector must exist

    const targetId = mouseWheelControlTargetSelect.value;
    const slider = document.getElementById(targetId);

    if (!slider || slider.type !== 'range') {
        return; // Allow page scrolling if the target isn't a valid range slider
    }

    // Check if target slider or its control group is visible
    const sliderGroup = slider.closest('.control-group');
    const style = window.getComputedStyle(sliderGroup || slider);
    if (!style || style.display === 'none' || style.visibility === 'hidden') {
        return; // Allow page scrolling if target is hidden
    }

    // Prevent page scrolling ONLY when wheeling over a valid, visible slider target
    event.preventDefault();
    event.stopPropagation();

    // Calculation Logic
    const delta = event.deltaY;
    const direction = delta < 0 ? 1 : -1;

    const step = parseFloat(slider.step) || 0.01;
    const min = parseFloat(slider.min);
    const max = parseFloat(slider.max);
    const currentValue = parseFloat(slider.value);
    const range = max - min;

    // Adjust step size based on parameter sensitivity or slider range
    let effectiveStep = step;
    if (range > 500 && step < 10) { // Larger steps for wide ranges like frequency/cutoff
       effectiveStep = Math.max(step * 5, range / 100);
    } else if (step < 0.01) { effectiveStep *= 3;
    } else if (step < 0.1) { effectiveStep *= 2; }

    // Increase sensitivity with Shift key
    if (event.shiftKey) { effectiveStep *= 0.2; } // Finer control with Shift

    let newValue = currentValue + direction * effectiveStep;

    // Clamp the value within min/max
    newValue = Math.max(min, Math.min(max, newValue));

    // Round to slider's step precision
    const numDecimalPlaces = (step.toString().split('.')[1] || '').length;
    if (numDecimalPlaces > 0) {
        newValue = parseFloat(newValue.toFixed(numDecimalPlaces));
    } else {
        newValue = Math.round(newValue);
    }

    // Apply the change if it's different
    if (newValue !== currentValue) {
        slider.value = newValue;
        slider.dispatchEvent(new Event('input', { bubbles: true }));
    }
}

// --- Initial Setup Function (Core UI/Audio, excluding MIDI init) ---
function initializeControlsAndAudio() {
    // Get canvas contexts & setup canvases
    // References are global, get contexts here
    if (oscilloscopeCanvas) oscilloscopeCtx = oscilloscopeCanvas.getContext('2d');
    if (spectrumCanvas) spectrumCtx = spectrumCanvas.getContext('2d');

    function setupCanvas(canvas, ctx) {
        if (!canvas || !ctx) return;
        try {
            const dpr = window.devicePixelRatio || 1;
            const rect = canvas.getBoundingClientRect();
            // Check for valid dimensions before setting
            if (rect.width > 0 && rect.height > 0) {
                canvas.width = rect.width * dpr;
                canvas.height = rect.height * dpr;
                ctx.scale(dpr, dpr); // Scale context for High DPI
                canvas.style.width = `${rect.width}px`;
                canvas.style.height = `${rect.height}px`;
            } else {
                // console.warn(`Canvas ${canvas.id} has zero dimensions.`);
            }
        } catch (e) {
            console.error(`Error setting up canvas ${canvas?.id}:`, e);
        }
    }
    setupCanvas(oscilloscopeCanvas, oscilloscopeCtx);
    setupCanvas(spectrumCanvas, spectrumCtx);

    // Initialize UI text values (Check elements exist)
    if(freqValueSpan && frequencySlider) freqValueSpan.textContent = frequencySlider.value;
    if(pluckinessValueSpan && pluckinessSlider) pluckinessValueSpan.textContent = pluckinessSlider.value;
    if(modalPresetASelect) modalPresetASelect.value = 'bell';
    if(modalPresetBSelect) modalPresetBSelect.value = 'levitationChamber';
    if(modalMorphSlider && modalMorphValueSpan) { modalMorphSlider.value = 0.0; modalMorphValueSpan.textContent = parseFloat(modalMorphSlider.value).toFixed(2); }
    if(chaosFeedbackSlider && chaosFeedbackValueSpan) { chaosFeedbackSlider.value = 0.0; chaosFeedbackValueSpan.textContent = parseFloat(chaosFeedbackSlider.value).toFixed(2); }
    if(delayTimeValueSpan && delayTimeSlider) delayTimeValueSpan.textContent = delayTimeSlider.value;
    if(delayFeedbackValueSpan && delayFeedbackSlider) delayFeedbackValueSpan.textContent = delayFeedbackSlider.value;
    if(filterCutoffValueSpan && filterCutoffSlider) filterCutoffValueSpan.textContent = filterCutoffSlider.value;
    if(masterGainValueSpan && masterGainSlider) masterGainValueSpan.textContent = masterGainSlider.value;
    if(mouseWheelControlTargetSelect) mouseWheelControlTargetSelect.value = 'frequency';
    if(panLfoRateValueSpan && panLfoRateSlider) panLfoRateValueSpan.textContent = parseFloat(panLfoRateSlider.value).toFixed(2);
    if(panLfoDepthValueSpan && panLfoDepthSlider) panLfoDepthValueSpan.textContent = parseFloat(panLfoDepthSlider.value).toFixed(2);
    if(grainAmountValueSpan && grainAmountSlider) grainAmountValueSpan.textContent = parseFloat(grainAmountSlider.value).toFixed(3);

    // Set initial engine & update UI
    currentEngine = 'string';
    if(engineStringRadio) engineStringRadio.checked = true;
    updateUIForEngine(); // Sets labels/ranges/values + hides/shows controls

    // Set initial Audio Param values (Check sliders exist)
    const now = audioContext.currentTime;
    if(filterCutoffSlider) masterFilter.frequency.setValueAtTime(parseFloat(filterCutoffSlider.value), now);
    if(delayTimeSlider) delay.delayTime.setValueAtTime(parseFloat(delayTimeSlider.value), now); // Use setValueAtTime for consistency
    if(delayFeedbackSlider) delayFeedback.gain.setValueAtTime(parseFloat(delayFeedbackSlider.value), now);
    if(masterGainSlider) masterVolume.gain.setValueAtTime(parseFloat(masterGainSlider.value), now);
    // LFO initial params already set where LFOs are created

    // Assign engine-specific classes to parent control groups
    document.getElementById('pluckiness')?.closest('.control-group')?.classList.add('engine-specific', 'string-specific', 'hybrid-specific');
    modalSpecificDiv?.classList.add('engine-specific', 'modal-specific', 'hybrid-specific');

     // Start the visualization loop
     if (oscilloscopeCtx || spectrumCtx) {
        stopVisualize(); // Stop previous loop if any
        visualize();
     } else {
        console.warn("Visualizations disabled (canvas contexts not found).");
     }

     // --- MIDI Initialization is MOVED outside this function, into the DOMContentLoaded listener ---
}


// --- WAIT FOR DOM BEFORE INITIALIZING AND SETTING UP LISTENERS/MIDI ---
document.addEventListener('DOMContentLoaded', () => {
    console.log("DOM fully loaded and parsed");

    // --- Initialize Core Synth UI & Audio ---
    initializeControlsAndAudio(); // Run the main init function

    // --- Setup Event Listeners for UI Controls ---
    setupEventListeners(); // Attach listeners to sliders, buttons, etc.

    // --- Setup MIDI AFTER DOM is ready and core init is done ---
    // Use the global references to MIDI DOM elements (already fetched at top)
    if (navigator.requestMIDIAccess) {
        console.log('Attempting to access MIDI...');
        navigator.requestMIDIAccess({ sysex: false }) // sysex false is safer
            .then(onMIDISuccess, onMIDIFailure); // Callbacks handle UI updates
    } else {
        console.warn('WebMIDI API not supported in this browser.');
        if (midiStatus) { // Check if element exists before using it
            midiStatus.textContent = 'MIDI not supported by your browser.';
            midiStatus.style.color = 'var(--accent-color-warning)';
        }
        if (midiDeviceSelector) { // Check if element exists
             midiDeviceSelector.disabled = true;
        }
    }

    // Add listener for MIDI device selection change ONLY if the element exists
    if (midiDeviceSelector) {
        midiDeviceSelector.addEventListener('change', (event) => {
            connectMIDIInput(event.target.value);
        });
         console.log("MIDI selector 'change' listener added.");
    } else {
         // This error suggests a problem with the HTML or element ID
         console.error("CRITICAL: Could not find MIDI device selector element (id='midiInDeviceSelector')!");
         if (midiStatus) {
              midiStatus.textContent = 'Error: MIDI UI Element missing.';
              midiStatus.style.color = 'var(--accent-color-warning)';
         }
    }

    console.log("Advanced Synth Initialization Complete.");
    if (audioContext.state === 'suspended') {
        console.log("AudioContext is suspended. Requires user interaction (click/key) to start audio.");
        // Optionally display a message to the user on the page
        // e.g., document.getElementById('audio-status-message').textContent = "Click anywhere to enable audio";
    } else {
         console.log("AudioContext state:", audioContext.state);
    }
});


// --- Computer Keyboard Event Listeners (Attach directly to window) ---
window.addEventListener('keydown', async (event) => {
    // Ignore if modifier keys are pressed or if it's a repeat event
    if (event.metaKey || event.ctrlKey || event.altKey || event.repeat) return;
    // Don't trigger computer keys if typing in an input/select/textarea
    const activeTag = document.activeElement?.tagName;
    if (activeTag === 'INPUT' || activeTag === 'SELECT' || activeTag === 'TEXTAREA') return;

    const frequency = getFrequencyFromKey(event.key);
    if (frequency !== null) { // Check if the key is mapped
        event.preventDefault(); // Prevent default actions like typing or page scroll
        if (audioContext.state === 'suspended') await resumeAudioContext();

        // Only play if audio context resumed successfully
        if(isAudioContextResumed) {
            const keyLower = event.key.toLowerCase();
            // Check pressedKeys *before* playing to handle retriggering correctly via playNote
            if (!pressedKeys.has(keyLower)) {
                // playNote now handles adding to pressedKeys
                playNote(keyLower, frequency, 90); // Use character key for tracking, default velocity
            } else {
                 // If key is already in pressedKeys, call playNote anyway to handle retrigger logic
                 playNote(keyLower, frequency, 90);
            }
        }
    }
});

window.addEventListener('keyup', (event) => {
    const keyLower = event.key.toLowerCase();
    if (keyToNoteMap.hasOwnProperty(keyLower)) {
         // Stop note using the character as identifier
         // stopNote now handles removing from pressedKeys
         stopNote(keyLower);
    }
});