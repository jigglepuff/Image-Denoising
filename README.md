# Image-Denosing
Summary: Sample codes to denoise a noisy image

Input image: .RAW format

Programs: 
1) basic_denoise_filters.cpp + Pixel3.h - Implementation of i) Linear/Averaging filter, ii) Non-linear/ median filter, 
                                          iii) Gaussian filter on RGB intensities, iv) Bilateral Gaussian filter on spatial position and RGB intensities
2) non_local_mean.cpp + Pixel3.h - Inspired by an academic paper (see http://www.coe.utah.edu/~cs7640/readings/NL-means.pdf). 
                                   The idea is to take guassian function of the neighborhood's local average values to produce better denosing results.  
