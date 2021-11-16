# Thesis
Performance analysis of LoRa networks with Successive Interference Cancellation.

# Abstract
The increase of the connected Internet-of-Things (IoT) devices, necessitates the development
of low-power wide-area networks (LPWANs). LoRa technology is a popular and promising solution
for communication in LPWANs. This thesis studies the use of Successive Interference Cancellation
(SIC) to decode received signals in LoRa networks. The proposed method utilises the capture effect
and designs a new receiver structure, capable of retrieving collided packets using SIC. A complete
theoretical analysis is provided and closed-form expressions are derived for the successful 
decoding of packets via SIC taking path loss, fading, noise, and interference into account. The
theoretical results are validated by Monte Carlo simulations.

# Problem Description
A circular coverage region is considered (R km radius) where there is only one gateway and multiple end devices. 
The baseline system's model performance is given in the figures below for two different cell radius, 6 km and 12 km.
<p float="left">
  <img src="/plots/basic_model.png" width="400" />
  <img src="/plots/basic_model2.png" width="400" /> 
</p>

In the figures solid lines are obtained via numerical evaluation and, whereas the markers represent simulation results.
- Connection Probability: A node is connected to the gateway if the SNR (Signal-to-Noise Ratio) is above a default threshold.
- Capture Probability: For a successful reception in the presence of interference, the SIR (Signal-to-Interference Ratio) must exceeds the capture threshold.
- Coverage Probability: A node is in coverage if it is connected to the gateway and there is no collision.

# Results

SIC technique is used to decode only one interfering message, even if there are more and system's performance is evaluated under specific circumstances. Since SIC is a method for mitigating interfernce, the connection probability will not be affected. 

## Capture probability
Applying SIC increases the capture probability up to 23.6%.
<p float="left">
  <img src="/plots/capture_prob.png" width="400" />
  <img src="/plots/capture_prob2.png" width="400" /> 
</p>

## Coverage probability
SIC increases the coverage probability up to 20% when the system radius is 6 km. A smaller improvement is observed at 12 km, since the noise (connection probability) is the primary cause of drop in performance.
<p float="left">
  <img src="/plots/coverage.png" width="400" />
  <img src="/plots/coverage2.png" width="400" /> 
</p>

# How to run

# Requirments

# License
