# eDNA-Latent-Copies
Models for eDNA concentration estimation that treat copy number as latent count variables observed with error (nondetection and continuous measurement error given detection)

Basic files to simulate from and fit models that regard the number of copies in each PCR rep as Poisson random variables. Other count distributions could be used (e.g., negative binomial), if appropriate.

Most files are for an eDNA release experiment. Minimal info for now, see files for the abstract and model diagrams.

See "Temporal AR1 testscript.R" for a different process model describing the concentration at 1 site through time with temporal correlation. Observation model is the same as for the release experiment models.

Inhibitor models may take a long time to converge, or not converge at all if the data is not sufficient for the model (leaving this vague for now).