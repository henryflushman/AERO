import numpy as np
import statistics as stat
import math
import scipy.stats as stats

t = [75.4, 74.1, 73.7, 74.1, 74.7, 74.8, 75.8, 73.8, 76.0, 74.1, 74.0, 74.2, 75.0, 75.8, 74.8, 73.7, 74.3, 74.5, 74.2, 75.2, 75.5]

mean = sum(t)/len(t)
t_sort = sorted(t)
median = t_sort[round(len(t)/2)]
counts = {}
for num in t:
    counts[num] = counts.get(num, 0) + 1

mode = max(counts, key=counts.get)

n = len(t)
ss = sum((ti - mean) ** 2 for ti in t)
variance = ss / (n - 1)

std_dev = math.sqrt(variance)

alpha = 0.05

df = n - 1
chi2_lower = stats.chi2.ppf(alpha / 2, df)
chi2_upper = stats.chi2.ppf(1 - alpha / 2, df)

a = math.sqrt((df * std_dev**2) / chi2_upper)
b = math.sqrt((df * std_dev**2) / chi2_lower)

min_temp = -20
max_temp = 100
measured_temp = 22.3
uncertainty_percent = 8
full_range = max_temp - min_temp
uncertainty = (uncertainty_percent / 100) * full_range
coldest_temp = measured_temp - uncertainty

print(f"n = {n}")
print(f"mean = {mean:.4f}")
print(f"median = {median:.4f}")
print(f"mode = {mode}")
print(f"sample variance = {variance:.6f}")
print(f"sample std dev = {std_dev:.6f}")
print(f"95% CI for sigma = [{a:.6f}, {b:.6f}]")
print(f"coldest_temp (with {uncertainty_percent}% of full-range uncertainty) = {coldest_temp:.2f}")