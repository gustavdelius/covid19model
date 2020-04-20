# Bayesian inference of the COVID19 Infection Fatality Ratio from mortality data

We use the model by the Imperial group at <https://github.com/ImperialCollegeLondon/covid19model>
to infer the ratio between the number of infections and the number of deaths for
COVID-19. This work is described in the report at 
[`covid19_IFR_report.pdf`](https://github.com/gustavdelius/covid19model/blob/master/covid19_IFR_report.pdf).
Here is the abstract:

**Abstract:** We use an established semi-mechanistic Bayesian hierarchical model of the COVID-19 pandemic, driven by European mortality data, to 
 estimate the prevalence of immunity. We allow the infection-fatality ratio (IFR) to vary, adapt the model's priors to better reflect emerging information, and  re-evaluate the model fitting in the light of current mortality data. The results indicate that the IFR of COVID-19 may be an order of magnitude smaller than the current consensus, with the corollary that the virus is more prevalent than currently believed. These results emerge from a simple model and ought to be treated with caution. They emphasise the value of rapid community-scale antibody testing when this becomes available.

You can find the Stan model to reproduce the results from that report in the
file [`stan-models/step_ifr.stan`](https://github.com/gustavdelius/covid19model/blob/master/stan-models/step_ifr.stan).

Figures with the results from the data up to 18/04/2020 for 14 European countries 
can be found in [`figures/18_04_step_ifr`](https://github.com/gustavdelius/covid19model/tree/master/figures/18_04_step_ifr).