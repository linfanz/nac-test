# nac-test

Simulations for the NAC family of tests. See the [paper](https://arxiv.org/abs/2012.15047) for more details.

The code relies on the [`nett`](https://github.com/aaamini/nett) package.

## simulations 

```model selection.R```: Compare different methods' accuracy in selecting the correct community numbers of DCSBM. Including 3 different designs of the connectivity matrix.

```ROC.R```: Generate ROC curves

## FB-100 analysis

FB-100 data has not been uploaded yet. ```FBhist``` contains histogram profiles of the networks. ```FBoutput``` contains other output of the networks.

```FB_hist.R```: create a histogram profile of SNAC+ statistics applied to FB networks and DCSBM with planted partition structure.

```FB_hist_reduced.R```: create a histogram profile of SNAC+ statistics applied to FB sub networks (with degree larger than 75 percentile) and DCSBM with planted partition structure.

```FB_norm_stat_plot.R```: normalized plots of several statistics applied to networks

```FB_smooth_profile.R```: create smooth curve of repeated SNAC+ statistics on networks

```FB_net_plot.R```: visualize networks with nodes colored using labels


