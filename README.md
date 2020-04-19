# Bayesian inference of the COVID19 Infection Fatality Ratio from mortality data

We use the model by the Imperial group at <https://github.com/ImperialCollegeLondon/covid19model>
to infer the ratio between the number of infections and the number of deaths for
COVID-19. This work is described in the report at ...

You can find the Stan model to reproduce the results from that report in the
file `stan-models/step_ifr.stan`.

Figures with the results from the data up to 18/04/2020 for 14 European countries 
can be found in `figurs/18_04_step_ifr`.