# Show energy compaction to motivate vector CS
# Must be done later: normalize the different coefficients to have same sum(?)
# plot all three on a chart with logarithmic x axis.

# See https://matplotlib.org/examples/pylab_examples/log_demo.html

# Lena is from http://links.uwaterloo.ca/Repository.html

import matplotlib.pyplot as plt

from scipy.misc import imread, imsave
from scipy.fftpack import dct, idct
import numpy as np
import pywt

# 2D DCT from https://stackoverflow.com/questions/34890585
def dct2 (block):
	return dct(dct(block.T, norm = 'ortho').T, norm = 'ortho')

def idct2(block):
	return idct(idct(block.T, norm = 'ortho').T, norm = 'ortho')

lena = imread("lena2.tif", "L")
lena = lena/255.0
dct_lena = dct2(lena)

#print lena
#print a

q = dct_lena.flatten()
num_coeffs = len(q)

wavelet_type = "haar"
wave = pywt.wavedec2(lena, wavelet_type)
wave_flattened = []
# coerce the wavelet coefficients to one long list
# (this loses some structure)
for i in wave:
	wave_flattened.extend(list(np.array(i).flatten()))

sorted_plain = -np.sort(-np.abs(lena.flatten()))
sorted_dct_coeffs = -np.sort(-np.abs(q))
sorted_dwt_coeffs =  -np.sort(-np.abs(wave_flattened))

dct_sample_count = np.array(range(num_coeffs))+1
dwt_sample_count = np.array(range(len(sorted_dwt_coeffs)))+1

#plt.semilogx(x, sorted_dct_coeffs)

plt.subplot(211)
plt.title('Distribution of coefficients, Lena')
plt.plot(dct_sample_count, sorted_dct_coeffs, label="DCT-II")
plt.plot(dwt_sample_count, sorted_dwt_coeffs, label="Haar wavelet")
plt.legend()
plt.grid(True)

plt.subplot(212)
plt.title('Distribution of coefficients, Lena, log plot')
plt.semilogx(dct_sample_count, sorted_dct_coeffs, label="DCT-II")
plt.semilogx(dwt_sample_count, sorted_dwt_coeffs, label="Haar wavelet")
plt.grid(True)
plt.legend()
plt.show()

# Do a simple 1:50 sparsification in both domains, just discarding all 
# values below this threshold

dct_threshold = abs(sorted_dct_coeffs[len(sorted_dct_coeffs)/50])
print dct_threshold
thresholded_w_dct = np.copy(dct_lena)
thresholded_w_dct = thresholded_w_dct * (abs(thresholded_w_dct) > dct_threshold)
lena_dct_recoded = idct2(thresholded_w_dct) * 255.0
imsave("lenadct.tif", lena_dct_recoded)

#pywt.threshold(thresh_dct, dct_threshold, mode="soft", substitute=dct_threshold)

dwt_threshold = abs(sorted_dwt_coeffs[len(sorted_dwt_coeffs)/50])
thresholded = [pywt.threshold(i, dwt_threshold, "hard")  for i in wave]

lena_dwt_recoded = pywt.waverec2(thresholded, wavelet_type) * 255.0
print lena_dwt_recoded
imsave("lenadwt.tif", lena_dwt_recoded)
