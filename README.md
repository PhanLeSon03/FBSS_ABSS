Dependent: install RIR-Generator (https://github.com/ehabets/RIR-Generator)

`Beam_coeff7mics.m`: generated beamformer coefficients 

`FBSS.m`: fixed beamforming with suppressed sidelobe

`Beam_SS.m`: waveform comparisons, free field, two inferences

`Beam_SS_RIR.m`: waveform comparisons including RIR and DF, set `diffuse_noise` to 0/1 for select white noise or diffuse noise

`Beam_SS_SIR.m`: beamforming with suppresed sidelobe, oSINR over SIR evaluations between different approaches:GSC, FBSS at SNR = 0 and -10dB

`Beam_SS_SIR_RIR.m`: beamforming with suppresed sidelobe, SIR evaluations between different approaches, room impulse response:GSC, ABSS, FBSS at SNR = 10 and 20dB
