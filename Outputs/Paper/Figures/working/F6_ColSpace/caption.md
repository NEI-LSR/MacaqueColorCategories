# colorspace_everySecond
Colors uniformally sampled in CIELUV, plotted in CIELUV

# combined_TCC-0att_fullremap-workspace_230510_behaviorally-derived-colorspace_everySecond
Colors uniformally sampled in CIELUV, spatial relationships between colors determined from the similarity structure uncovered from the behavioral data. 
Though the sampling interval is nominally unifrorm, in reality, the colors in the blue and red areas are more similar to their neighbours than the colors in the green and lavender areas are to their neighbours.

# newEqualSampling
Colors uniformally sampled in the new behaviorally derived colorspace, plotted in CIELUV. 
You can see that if we were using CIELUV to sample, we should really sample less points from blue/red if we want a "uniform" sampling.

# colspace_ssnuModelOfSGData_behaviorally-derived-colorspace_everySecond (and sm_ssnuModelOfSGData)
Behavioral data was simulated using a hypothetical perfectly uniform colorspace, but with cognitive biases. 
Then a ssnu model was fit to the data, and a colorspace recovered.
This shows that by fitting the wrong model to the data, it can look like there is a non-uniform colorspace, when there is actually not.
