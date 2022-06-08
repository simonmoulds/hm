Overview
========

Environmental models need a large amount of boilerplate code to handle aspects such as I/O, model configuration, and timestepping. Model developers tend to implement their own bespoke framework to perform these tasks. This results in excessive duplication across environmental models, and makes it harder to implement multi-model experiments. The goal of **hm** is twofold: to simplify the model development process and allow scientists to focus on the science code, and to provide a common interface to environmental models which makes it easier to implement multi-model ensembles and model coupling.
