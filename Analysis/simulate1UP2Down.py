
# %% 
import numpy as np
from collections import Counter
from itertools import product

import matplotlib.pyplot as plt

def simulate_1up2down(n_trials=60, p_correct=0.71, start_amp=1.0, step_size=0.1, min_amp=0.1, max_amp=2.0):
    """
    Simulate a 1-up 2-down staircase procedure.
    Args:
        n_trials: Number of trials to simulate.
        p_correct: Target probability of correct response (converges to ~0.71 for 1up2down).
        start_amp: Starting amplitude.
        step_size: Step size for amplitude adjustment.
        min_amp: Minimum amplitude.
        max_amp: Maximum amplitude.
    Returns:
        amps: List of amplitude values per trial.
        responses: List of simulated responses (1=correct, 0=incorrect).
    """
    amps = [start_amp]
    responses = []
    correct_count = 0

    for i in range(n_trials):
        # Simulate response: 1=correct, 0=incorrect
        resp = np.random.rand() < p_correct
        responses.append(int(resp))

        # 1-up 2-down logic
        if resp:
            correct_count += 1
            if correct_count == 2:
                # Decrease amplitude after 2 consecutive correct
                new_amp = max(min_amp, amps[-1] - step_size)
                amps.append(new_amp)
                correct_count = 0
            else:
                amps.append(amps[-1])
        else:
            # Increase amplitude after 1 incorrect
            new_amp = min(max_amp, amps[-1] + step_size)
            amps.append(new_amp)
            correct_count = 0

    return amps[:-1], responses

# Parameters
n_trials = 60
p_correct = 0.71  # 1up2down converges to ~0.71
start_amp = 1.0
step_size = 0.1

# Simulate
amps, responses = simulate_1up2down(n_trials, p_correct, start_amp, step_size)

# Plot
plt.figure(figsize=(10, 5))
plt.plot(amps, marker='o', label='Amplitude')
plt.xlabel('Trial')
plt.ylabel('Amplitude')
plt.title('1-up 2-down Staircase Adjustment Curve')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
# %%
from itertools import product 
N= 20
# 生成所有可能的10次反应组合（0=错误, 1=正确）
all_combos = list(product([0, 1], repeat=N))
# 计算每种正确数出现的概率
prob_dict = {}
p = 0.71
for combo in all_combos:
    correct = sum(combo)
    prob = (p ** correct) * ((1 - p) ** (N - correct))
    prob_dict[correct] = prob_dict.get(correct, 0) + prob

for correct in range(N + 1):
    cum_prob = sum(prob_dict[i] for i in range(correct + 1))
    print(f"{correct} correct: probability = {prob_dict[correct]:.6f}, cumulative probability = {cum_prob:.6f}")
    
# %%
from scipy import stats
# Calculate the probability of getting between 66% and 76% correct in 12 trials
stats.beta.cdf(0.76, 281, 121) - stats.beta.cdf(0.66, 281, 121)

# %%
