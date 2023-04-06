#!/usr/bin/env python
# coding: utf-8

# # MadMiner physics tutorial (part 2A)
# 
# Johann Brehmer, Felix Kling, Irina Espejo, and Kyle Cranmer 2018-2019

# In this second part of the tutorial, we'll generate events and extract the observables and weights from them. You have two options: In this notebook we'll do this at parton level, in the alternative part 2b we use Delphes.

# ## 0. Preparations

# Before you execute this notebook, make sure you have a running installation of MadGraph.

# In[1]:


import os
import logging
import numpy as np

from madminer.core import MadMiner
from madminer.lhe import LHEReader
from madminer.sampling import combine_and_shuffle
from madminer.plotting import plot_distributions
from particle import Particle


# In[2]:


# MadMiner output
logging.basicConfig(
    format="%(asctime)-5.5s %(name)-20.20s %(levelname)-7.7s %(message)s",
    datefmt="%H:%M",
    level=logging.DEBUG,
)

# Output of all other modules (e.g. matplotlib)
for key in logging.Logger.manager.loggerDict:
    if "madminer" not in key:
        logging.getLogger(key).setLevel(logging.WARNING)


# In[3]:


#mg_dir = os.getenv("/Users/zimermma/work/MG5_aMC_v3_1_0/")
mg_dir="/nfs/dust/cms/user/zimermma/madminer/MG5_aMC_v2_8_3_2"

# ## 1. Generate events

# Please enter here the environment variable pointing to your MG5 installation folder.

# Let's load our setup:

# In[7]:


miner = MadMiner()
miner.load("data/setup.h5")


# In a next step, MadMiner starts MadGraph to generate events and calculate the weights. You can use `run()` or `run_multiple()`; the latter allows to generate different runs with different run cards and optimizing the phase space for different benchmark points. 
# 
# In either case, you have to provide paths to the process card, run card, param card (the entries corresponding to the parameters of interest will be automatically adapted), and an empty reweight card. Log files in the `log_directory` folder collect the MadGraph output and are important for debugging.
# 
# The `sample_benchmark` (or in the case of `run_all`, `sample_benchmarks`) option can be used to specify which benchmark should be used for sampling, i.e. for which benchmark point the phase space is optimized. If you just use one benchmark, reweighting to far-away points in parameter space can lead to large event weights and thus large statistical fluctuations. It is therefore often a good idea to combine a lot of events at the "reference hypothesis" (for us the SM) and smaller samples from other benchmarks that span the parameter space.
# 

# In[8]:


miner.run(
    sample_benchmark="sm",
    mg_directory=mg_dir,
    mg_process_directory="./mg_processes/signal1",
    proc_card_file="cards/proc_card_signal.dat",
    param_card_template_file="cards/param_card_template.dat",
    run_card_file="cards/run_card_signal_large.dat",
    log_directory="logs/signal",
    python_executable="python",
)


# In[9]:


additional_benchmarks = ["ctG", "neg_ctG", "ctu8", "neg_ctu8","morphing_basis_vector_5"]


# In[10]:


miner.run_multiple(
    sample_benchmarks=additional_benchmarks,
    mg_directory=mg_dir,
    mg_process_directory="./mg_processes/signal2",
    proc_card_file="cards/proc_card_signal.dat",
    param_card_template_file="cards/param_card_template.dat",
    run_card_files=["cards/run_card_signal_small.dat"],
    log_directory="logs/signal",
    python_executable="python",
)


# This will take a moment -- time for a coffee break!
# 
# After running any event generation through MadMiner, you should check whether the run succeeded: are the usual output files there, do the log files show any error messages? MadMiner does not (yet) perform any explicit checks, and if something went wrong in the event generation, it will only notice later when trying to load the event files.

# ### Backgrounds

# We can also easily add other processes like backgrounds. An important option is the `is_background` keyword, which should be used for processes that do *not* depend on the parameters theta. `is_background=True` will disable the reweighting and re-use the same weights for all cross sections.
# 
# To reduce the runtime of the notebook, the background part is commented out here. Feel free to activate it and let it run during a lunch break.

# In[26]:


#miner.run(
#    is_background=True,
#    sample_benchmark='sm',
#    mg_directory=mg_dir,
#    mg_process_directory='./mg_processes/background',
#    proc_card_file='cards/proc_card_background.dat',
#    param_card_template_file='cards/param_card_template.dat',
#    run_card_file='cards/run_card_background.dat',
#    log_directory='logs/background',
#)


# Finally, note that both `MadMiner.run()` and `MadMiner.run_multiple()` have a `only_create_script` keyword. If that is set to True, MadMiner will not start the event generation directly, but prepare folders with all the right settings and ready-to-run bash scripts. This might make it much easier to generate Events on a high-performance computing system. 

# In[34]:





# ## 2. Prepare analysis of the LHE samples

# The `madminer.lhe` submodule allows us to extract observables directly from the parton-level LHE samples, including an approximate description of the detector response with smearing functions. The central object is an instance of the `LHEProcessor` class, which has to be initialized with a MadMiner file:

# In[11]:


lhe = LHEReader("data/setup.h5")


# After creating the `LHEReader` object, one can add a number of event samples (the output of running MadGraph in step 1) with the `add_sample()` function.
# 
# In addition, you have to provide the information which sample was generated from which benchmark with the `sampled_from_benchmark` keyword, and set `is_background=True` for all background samples.

# In[12]:


lhe.add_sample(
    lhe_filename="mg_processes/signal1/Events/run_01/unweighted_events.lhe.gz",
    sampled_from_benchmark="sm",
    is_background=False,
    k_factor=1.0,
)
for i, benchmark in enumerate(additional_benchmarks):
    lhe.add_sample(
        lhe_filename="mg_processes/signal2/Events/run_0{}/unweighted_events.lhe.gz".format(i + 1),
        sampled_from_benchmark=benchmark,
        is_background=False,
        k_factor=1.0,
    )


#lhe.add_sample(
#    lhe_filename='mg_processes/background/Events/run_01/unweighted_events.lhe.gz',
#    sampled_from_benchmark='sm',
#    is_background=True,
#    k_factor=1.0,
#)


# ## 3. Smearing functions to model the detector response

# Now we have to define the smearing functions that are used (in lieu of a proper shower and detector simulation). Here we will assume a simple 10% uncertainty on the jet energy measurements and a $\pm 0.1$ smearing for jet $\eta$ and $\phi$. The transverse momenta of the jets are then derived from the smeared energy and the on-shell condition for the quarks (this is what `pt_resolution_abs=None` does). The photons from the Higgs are assumed to be measured perfectly (otherwise we'd have to call `set_smearing` another time with `pdgis=[22]`).

# In[13]:


# Partons giving rise to jets
particles = [
    *Particle.findall(lambda p: p.pdgid.is_quark),
    *Particle.findall(pdg_name="g"),
]

lhe.set_smearing(
    pdgids=[int(p.pdgid) for p in particles],
    energy_resolution_abs=0.0,
    energy_resolution_rel=0.1,
    pt_resolution_abs=None,
    pt_resolution_rel=None,
    eta_resolution_abs=0.1,
    eta_resolution_rel=0.0,
    phi_resolution_abs=0.1,
    phi_resolution_rel=0.0,
)


# In addition, we can define noise that only affects MET. This adds Gaussian noise with mean 0 and std `abs_ + rel * HT` to MET_x and MET_y separately.

# In[14]:


lhe.set_met_noise(abs_=10.0, rel=0.05)


# ## 4. Observables and cuts

# The next step is the definition of observables, either through a Python function or an expression that can be evaluated. Here we demonstrate the latter, which is implemented in `add_observable()`. In the expression string, you can use the terms `j[i]`, `e[i]`, `mu[i]`, `a[i]`, `met`, where the indices `i` refer to a ordering by the transverse momentum. In addition, you can use `p[i]`, which denotes the `i`-th particle in the order given in the LHE sample (which is the order in which the final-state particles where defined in MadGraph).
# 
# All of these represent objects inheriting from scikit-hep [LorentzVectors](http://scikit-hep.org/api/math.html#vector-classes), see the link for a documentation of their properties. In addition, they have `charge` and `pdg_id` properties.
# 
# `add_observable()` has an optional keyword `required`. If `required=True`, we will only keep events where the observable can be parsed, i.e. all involved particles have been detected. If `required=False`, un-parseable observables will be filled with the value of another keyword `default`.
# 
# In a realistic project, you would want to add a large number of observables that capture all information in your events. Here we will just define two observables, the transverse momentum of the leading (= higher-pT) jet, and the azimuthal angle between the two leading jets.

# In[15]:


lhe.add_observable(
    "pt_j1",
    "j[0].pt",
    required=False,
    default=0.0,
)
lhe.add_observable(
    "delta_phi_jj",
    "j[0].deltaphi(j[1]) * (-1.0 + 2.0 * float(j[0].eta > j[1].eta))",
    required=True,
)
lhe.add_observable(
    "met",
    "met.pt",
    required=True,
)


# We can also add cuts, again in parse-able strings. In addition to the objects discussed above, they can contain the observables:

# In[16]:


lhe.add_cut("(a[0] + a[1]).m > 122.0")
lhe.add_cut("(a[0] + a[1]).m < 128.0")
lhe.add_cut("pt_j1 > 20.0")


# ## 5. Run analysis and store processes events

# The function `analyse_samples` then calculates all observables from the LHE file(s) generated before, applies the smearing, and checks which events pass the cuts:

# In[17]:


lhe.analyse_samples()


# The values of the observables and the weights are then saved in the HDF5 file. It is possible to overwrite the same file, or to leave the original file intact and save all the data into a new file as follows:

# In[42]:


lhe.save("data/lhe_data.h5")


# ## 6. Plot distributions

# Let's see what our MC run produced:

# In[46]:


_ = plot_distributions(
    filename="data/lhe_data.h5",
    parameter_points=["sm", np.array([10.0, 0.0])],
    line_labels=["SM", "BSM"],
    uncertainties="none",
    n_bins=20,
    n_cols=3,
    normalize=True,
    sample_only_from_closest_benchmark=True,
)


# ## 7. Combine and shuffle different samples

# To reduce disk usage, you can generate several small event samples with the steps given above, and combine them now. Note that (for now) it is essential that all of them are generated with the same setup, including the same benchmark points / morphing basis!
# 
# This is generally good practice even if you use just one sample, since the events might have some inherent ordering (e.g. from sampling from different hypotheses). Later when we split the events into a training and test fraction, such an ordering could cause problems.

# In[45]:


combine_and_shuffle(["data/lhe_data.h5"], "data/lhe_data_shuffled.h5")


# In[ ]:



MC_v2_8_3_2
