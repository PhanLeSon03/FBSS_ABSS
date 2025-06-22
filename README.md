## Dependencies

- [RIR-Generator](https://github.com/ehabets/RIR-Generator) — Required for simulating Room Impulse Responses (RIR).

## Script Descriptions

### Beamformer Design
- **`Beam_coeff7mics.m`**  
  Generates beamformer coefficients for a 7-microphone linear array.

### Core Experiments
- **`FBSS.m`**  
  Implements Fixed Beamforming with Suppressed Sidelobes (FBSS).

- **`Beam_SS.m`**  
  Adaptive Beamforming with Suppressed Sidelobes (ABSS) vs GSC:Compares beamformed waveforms under free-field conditions with two interference sources.


- **`Beam_SS_RIR.m`**  
  Compares waveforms using Room Impulse Responses (RIR) and Diffuse Field (DF) simulation.  
  - Use `diffuse_noise = 0` for white noise  
  - Use `diffuse_noise = 1` for diffuse noise

- **`Beam_SS_SIR.m`**  
  Evaluates output SINR (oSINR) over SIR across methods:  
  - GSC, ABSS, FBSS at SNR = 0 dB and –10 dB
  
- **`learning_rate.m`**  
  Evaluates output SINR (oSINR) over SIR across different learning rate values and approaches.

- **`Beam_SS_SIR_RIR.m`** *(Not used)*  
  oSINR/SIR evaluation under RIR conditions using GSC, ABSS, and FBSS at SNR = 10 dB and 20 dB.
