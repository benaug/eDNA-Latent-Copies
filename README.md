# eDNA-Latent-Copies
Models for eDNA concentration estimation that treat copy number as latent count variables observed with error (nondetection and continuous measurement error given detection)

Basic files to simulate from and fit models that regard the number of copies in each PCR rep as Poisson random variables. Other count distributions could be used (e.g., negative binomial), if appropriate.

These models are for an eDNA release experiment, but latent copy approach can be conditioned on other process models. Minimal info for now, see files for the abstract and model diagrams.

Inhibitor models may take a long time to converge, or not converge at all if the data is not sufficient for the model (leaving this vague for now).