# End-to-End 5G NR System

## Overview
This project implements an end-to-end 5G NR physical layer simulation system.
It covers the main signal processing blocks from transmitter to receiver, including modulation, MIMO-OFDM transmission, channel estimation, signal detection, and performance evaluation.
The goal is to study link-level behavior and analyze metrics such as BER, throughput, and capacity.

## Features
- End-to-end transmitter and receiver chain
- QAM modulation and demodulation
- OFDM processing
- MIMO transmission and detection
- Channel estimation and equalization
- LDPC encoding and decoding
- BER, throughput, and capacity analysis

## System Flow
Input image -> Bitstream -> Channel coding -> Modulation -> OFDM -> MIMO channel -> Channel estimation -> Detection -> Decoding -> Performance evaluation

## Project Structure
- `src/`: source code
- `results/`: simulation results
- `figures/`: generated plots
- `docs/`: project notes and references

## Requirements
- `MATLAB R2023a` or above
- Communication-related toolboxes if needed
- Custom `MEX` files if used in the project

## How to Run
1. Open the project in MATLAB.
2. Run the main simulation script.
3. Adjust parameters such as modulation order, SNR, MIMO size, and code rate.
4. Check output figures and result files in the corresponding folders.

## Notes
This repository is for research and academic study purposes.
