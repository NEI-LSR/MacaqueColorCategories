[intro]

In the task the monkey fixates, then is shown a cue colour, then that goes away, 4 choice colours come onto the screen, and then when the fixation disappears they're trained to select the choice which matches the colour of the cue for that trial.

They do this task very well, and we've collected data from 4 animals, about 200K trials in total.

A mixture modelling analysis recovers two attractor points, which can be interpreted as category centres: one blue/cool, and one orange/warm.

(Mixture model: a mixture between a uniform distribution (guessing) and a gaussian distribution centred on the correct answer. Here we add an additional parameter so that the gaussian doesn't need to be centred on the correct answer, and can instead show us whether the average response is biased.)

We realised however that there are two distinct potential sources for these biases.
The first being what people generally think of when they talk about bias: one colour is remembered as another.
The second can occur on this type of task if the distribution of stimuli is non-uniform with respect to the colorspace that is inherently used to complete the task. 
We sampled our stimuli from in the nominally uniform CIELUV colorspace, but there's no such thing as a colorspace that is truly uniform for every observer under every situation, and so we can assume that there *will* be at least some non-uniformity.
The knock on effect of this is that some colours will be more likely to be confused with neighbours on one side than on the other side (in the visual, the greenish-yellow at the center of the bullseye will be more likely to be confused with a green distractor than a yellow distractor). 
This would then be analysed as a bias towards green.

One way to distinguish these two types of bias is by using a similarity matrix: a matrix that lists the similarity between every colour, and every other colour.
(This can be thought of an extension to the TCC model's "similarity function", but rather than having a single similarity function that applies to all colours, we have separate similarity functions for each colour)
In the case of cognitive bias, we'd expect to see something like an s-curve around the diagonal - here, both greens and blues are pulled towards cyan. 
For the stimulus-space non-uniformity similarity matrix, based on the distorted colorspace illustrated above, we see an extension symmetrically around the diagonal in the area where the sampling is compressed (cyan), because colours in this area are closer to their neighbours than in other parts of the space. This will also result in biases towards the cyan.
Looking at data simulated using both of these very different similarity matrices we can see that they produce nearly identical results in a mixture model analysis.

If we fit generative models to our data (models that have access to the trial-by-trial data, and use different mechanisms to try to fit parameter values that best predict the data) we see that a model that allows for cognitive bias doesn't deliver an advantage over one without bias.
(Lower negative log likelihood = better)
We see however that a model that allows stimulus-space non-uniformity performs much better.
The fourth model that we fit has a single parameter for every single cell of the similarity matrix (we call it a free similarity matrix model), and this is the most flexible model - allowing for both cognitive bias *and* stimulus-space non-uniformity (and even different versions of both of the above). 
We see that it gets a lower NLL value (better), but is penalised for having so many parameters, and ends up with a higher AIC value.
Visually comparing the similarity matrices though, you can see that the stimulus-space non-uniformity matrix captures most of the large-scale structure of the free similarity matrix model.
Noteworthy as well is that the free similarity matrix appears to be broadly symmetrical around the diagonal, arguing against cognitive bias in this data.

From the stimulus-space non-uniformity model, we can extract the remapping of the stimuli. This gives us a representation of the latent behavioural colorspace. You can see compression around blue and around orange.

So far we've just been talking about the combined data across all 4 animals, but there are some individual differences.
M2 was noteworthy because his data showed evidence of an additional category in the yellow*, and if we look at a free similarity matrix model fit for his data alone (centred on that green category) we can see that it is driven by greens being pushed towards yellow, but curiously oranges aren't being attracted towards yellow, and so this appears to be a one-sided cognitive bias.

*(And technically a category in purple, but that's barely a zero crossing so we ignore that here)

[Conclusions]


