#!/usr/bin/env python
# coding: utf-8

# # MadMiner physics tutorial (part 3B)
#
# Johann Brehmer, Felix Kling, Irina Espejo, and Kyle Cranmer 2018-2019

# In part 3A of this tutorial we will finally train a neural network to estimate likelihood ratios. We assume that you have run part 1 and 2A of this tutorial. If, instead of 2A, you have run part 2B, you just have to load a different filename later.

# ## Preparations

# Make sure you've run the first tutorial before executing this notebook!

# In[1]:

import sys
import logging
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import argparse


#get_ipython().run_line_magic('matplotlib', 'inline')

from madminer.sampling import SampleAugmenter
from madminer import sampling
from madminer.ml import ScoreEstimator

parser = argparse.ArgumentParser(description="This is a script to calcualte the score of samples")

parser.add_argument(
    "-s", "--nsamples", type=int, default=500000,
    help="Total number of events to be drawn for training samples x ~ p(x|theta)"
    " as well as the joint score t(x, z|theta) "
)

parser.add_argument(
    "-e", "--nepochs", type=int, default=50,
    help="Number of Nepochs in the dnnn trianing "
)

parser.add_argument(
    "-l", "--hiddenlayers", type=int, default=2,
    help="Architecture of the hidden layers "
)

parser.add_argument(
    "-u", "--units", type=int, default=30,
    help="Number of layer units "
)

args = parser.parse_args()

DNNarch = (args.units,) * args.hiddenlayers

print("Starting the Score estimation with input parameters")
print("")
print("     Number of samples: ", args.nsamples)
print("     Number of epochs: ", args.nepochs)
print("     DNN architecture: ", DNNarch)
print("")

# In[2]:


# MadMiner output
logging.basicConfig(
    format="%(asctime)-5.5s %(name)-20.20s %(levelname)-7.7s %(message)s",
    datefmt="%H:%M",
    level=logging.INFO,
)

# Output of all other modules (e.g. matplotlib)
for key in logging.Logger.manager.loggerDict:
    if "madminer" not in key:
        logging.getLogger(key).setLevel(logging.WARNING)


# ## 1. Make (unweighted) training and test samples with augmented data

# At this point, we have all the information we need from the simulations. But the data is not quite ready to be used for machine learning. The `madminer.sampling` class `SampleAugmenter` will take care of the remaining book-keeping steps before we can train our estimators:
#
# First, it unweights the samples, i.e. for a given parameter vector `theta` (or a distribution `p(theta)`) it picks events `x` such that their distribution follows `p(x|theta)`. The selected samples will all come from the event file we have so far, but their frequency is changed -- some events will appear multiple times, some will disappear.
#
# Second, `SampleAugmenter` calculates all the augmented data ("gold") that is the key to our new inference methods. Depending on the specific technique, these are the joint likelihood ratio and / or the joint score. It saves all these pieces of information for the selected events in a set of numpy files that can easily be used in any machine learning framework.

# In[44]:


sampler = SampleAugmenter("data/lhe_data_Jonas_shuffled.h5")
# sampler = SampleAugmenter('data/delphes_data_shuffled.h5')


# The relevant `SampleAugmenter` function for local score estimators is `extract_samples_train_local()`. As in part 3a of the tutorial, for the argument `theta` you can use the helper functions `sampling.benchmark()`, `sampling.benchmarks()`, `sampling.morphing_point()`, `sampling.morphing_points()`, and `sampling.random_morphing_points()`.

# In[45]:


x, theta, t_xz, _ = sampler.sample_train_local(
    theta=sampling.benchmark("sm"),
    n_samples=args.nsamples,
    folder="./data/samples",
    filename="train_score",
)


# We can use the same data as in part 3a, so you only have to execute this if you haven't gone through tutorial 3a:

# In[5]:


_ = sampler.sample_test(
    theta=sampling.benchmark("sm"),
    n_samples=1000,
    folder="./data/samples",
    filename="test",
)


# ## 2. Train score estimator

# It's now time to build a neural network. Only this time, instead of the likelihood ratio itself, we will estimate the gradient of the log likelihood with respect to the theory parameters -- the score. To be precise, the output of the neural network is an estimate of the score at some reference parameter point, for instance the Standard Model. A neural network that estimates this "local" score can be used to calculate the Fisher information at that point. The estimated score can also be used as a machine learning version of Optimal Observables, and likelihoods can be estimated based on density estimation in the estimated score space. This method for likelihood ratio estimation is called SALLY, and there is a closely related version called SALLINO. Both are explained in ["Constraining Effective Field Theories With Machine Learning"](https://arxiv.org/abs/1805.00013) and ["A Guide to Constraining Effective Field Theories With Machine Learning"](https://arxiv.org/abs/1805.00020).
#
# The central object for this is the `madminer.ml.ScoreEstimator` class:

# In[50]:


estimator = ScoreEstimator(n_hidden=DNNarch)


# In[51]:


estimator.train(
    method="sally",
    x="data/samples/x_train_score.npy",
    t_xz="data/samples/t_xz_train_score.npy",
    n_epochs=args.nepochs,
)

estimator.save("models/sally")


# ## 3. Evaluate score estimator

# Let's evaluate the SM score on the test data

# In[31]:


estimator.load("models/sally")

t_hat = estimator.evaluate_score(x="data/samples/x_test.npy")


# Let's have a look at the estimated score and how it is related to the observables:

# In[60]:


x = np.load("data/samples/x_test.npy")

Nplots=10
fig = plt.figure(figsize=(12, 4*Nplots))

plt.rcParams.update({'font.size': 11, 'axes.labelsize': 13}) #change 14 to whatever font size you want

#mtt vs chel
ax = plt.subplot(Nplots, 2, 1)
sc = plt.scatter(
    x[:, 10],
    x[:, 0],
    c=t_hat[:, 0],
    s=25.0,
    cmap="BrBG",
    vmin=-0.8,
    vmax=+0.8,
)

cbar = plt.colorbar(sc)

cbar.set_label(r"$\hat{t}(x | \theta_{SM})$ for $c_{tG}$")
plt.xlabel(r"$c_{hel}$")
plt.ylabel(r"$m_{t\bart} [GeV]$")
plt.xlim(-1.05, 1.05)
plt.ylim(330, 1100)

ax = plt.subplot(Nplots, 2, 2)
sc = plt.scatter(
    x[:, 10],
    x[:, 0],
    c=t_hat[:, 1],
    s=25.0,
    cmap="BrBG",
    vmin=-0.65,
    vmax=+0.65,
)
cbar = plt.colorbar(sc)

cbar.set_label(r"$\hat{t}(x | \theta_{SM})$ for $c^8_{tq}$")
plt.xlabel(r"$c_{hel}$")
plt.ylabel(r"$m_{t\bart} [GeV]$")
plt.xlim(-1.05, 1.05)
plt.ylim(330, 1100)

#b2k vs chel
ax = plt.subplot(Nplots, 2, 3)
sc = plt.scatter(
    x[:, 5],
    x[:, 0],
    c=t_hat[:, 0],
    s=25.0,
    cmap="BrBG",
    vmin=-0.8,
    vmax=+0.8,
)

cbar = plt.colorbar(sc)

cbar.set_label(r"$\hat{t}(x | \theta_{SM})$ for $c_{tG}$")
plt.xlabel(r"$b^2_{k}$")
plt.ylabel(r"$m_{t\bart} [GeV]$")
plt.xlim(-1.05, 1.05)
plt.ylim(330, 1100)

ax = plt.subplot(Nplots, 2, 4)
sc = plt.scatter(
    x[:, 5],
    x[:, 0],
    c=t_hat[:, 1],
    s=25.0,
    cmap="BrBG",
    vmin=-0.65,
    vmax=+0.65,
)
cbar = plt.colorbar(sc)

cbar.set_label(r"$\hat{t}(x | \theta_{SM})$ for $c^8_{tq}$")
plt.xlabel(r"$b^2_{k}$")
plt.ylabel(r"$m_{t\bart} [GeV]$")
plt.xlim(-1.05, 1.05)
plt.ylim(330, 1100)

#thetaT vs chel
ax = plt.subplot(Nplots, 2, 5)
sc = plt.scatter(
    x[:, 3],
    x[:, 0],
    c=t_hat[:, 0],
    s=25.0,
    cmap="BrBG",
    vmin=-0.8,
    vmax=+0.8,
)

cbar = plt.colorbar(sc)

cbar.set_label(r"$\hat{t}(x | \theta_{SM})$ for $c_{tG}$")
plt.xlabel(r"$cos \Theta_t$")
plt.ylabel(r"$m_{t\bart} [GeV]$")
plt.xlim(-1.05, 1.05)
plt.ylim(330, 1100)

ax = plt.subplot(Nplots, 2, 6)
sc = plt.scatter(
    x[:, 3],
    x[:, 0],
    c=t_hat[:, 1],
    s=25.0,
    cmap="BrBG",
    vmin=-0.65,
    vmax=+0.65,
)
cbar = plt.colorbar(sc)

cbar.set_label(r"$\hat{t}(x | \theta_{SM})$ for $c^8_{tq}$")
plt.xlabel(r"$cos \Theta_t$")
plt.ylabel(r"$m_{t\bart} [GeV]$")
plt.xlim(-1.05, 1.05)
plt.ylim(330, 1100)

#mtt vs chel*
ax = plt.subplot(Nplots, 2, 7)
sc = plt.scatter(
    x[:, 11],
    x[:, 0],
    c=t_hat[:, 0],
    s=25.0,
    cmap="BrBG",
    vmin=-0.8,
    vmax=+0.8,
)

cbar = plt.colorbar(sc)

cbar.set_label(r"$\hat{t}(x | \theta_{SM})$ for $c_{tG}$")
plt.xlabel(r"$c_{hel*}$")
plt.ylabel(r"$m_{t\bart} [GeV]$")
plt.xlim(-1.05, 1.05)
plt.ylim(330, 1100)

ax = plt.subplot(Nplots, 2, 8)
sc = plt.scatter(
    x[:, 11],
    x[:, 0],
    c=t_hat[:, 1],
    s=25.0,
    cmap="BrBG",
    vmin=-0.65,
    vmax=+0.65,
)
cbar = plt.colorbar(sc)

cbar.set_label(r"$\hat{t}(x | \theta_{SM})$ for $c^8_{tq}$")
plt.xlabel(r"$c_{hel*}$")
plt.ylabel(r"$m_{t\bart} [GeV]$")
plt.xlim(-1.05, 1.05)
plt.ylim(330, 1100)

#yrap vs costheta**
ax = plt.subplot(Nplots, 2, 9)
sc = plt.scatter(
    x[:, 3],
    x[:, 1],
    c=t_hat[:, 0],
    s=25.0,
    cmap="BrBG",
    vmin=-0.8,
    vmax=+0.8,
)

cbar = plt.colorbar(sc)

cbar.set_label(r"$\hat{t}(x | \theta_{SM})$ for $c_{tG}$")
plt.xlabel(r"$cos \Theta_t$")
plt.ylabel(r"$y_{t}$")
plt.xlim(-1.05, 1.05)
plt.ylim(-4, 4)

ax = plt.subplot(Nplots, 2, 10)
sc = plt.scatter(
    x[:, 3],
    x[:, 1],
    c=t_hat[:, 1],
    s=25.0,
    cmap="BrBG",
    vmin=-0.65,
    vmax=+0.65,
)
cbar = plt.colorbar(sc)

cbar.set_label(r"$\hat{t}(x | \theta_{SM})$ for $c^8_{tq}$")
plt.xlabel(r"$cos \Theta_t$")
plt.ylabel(r"$y_{t}$")
plt.xlim(-1.05, 1.05)
plt.ylim(-4, 4)
plt.tight_layout()
#plt.show()


# In[76]:

figFormats = ['png','pdf','svg']
for figF in figFormats:
    fig.savefig("figures/3b_score_scatterplots_"+f"Nsamples_{str(args.nsamples)}_Nepochs_{str(args.nepochs)}_Nhlayers_{str(args.hiddenlayers)}_Units_{str(args.units)}."+figF,format=figF)


# In[ ]:
