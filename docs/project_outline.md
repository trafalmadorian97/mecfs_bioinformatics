

# Project Outline
## Goal
This project uses bioinformatic techniques to advance our understanding of ME/CFS and other poorly-understood diseases.


## Principles

When working on this project, there are a few key principles that govern how I write software and analyze data:



**Reproducibility**: It should be possible to reproduce any of main analyses with a few lines of Python code. 





**Verification**: As much as is feasible, code and analysis should be verified.

Verification techniques include:

- Implementing unit and integration tests to make sure that software continues to behave as expected, even as the software architecture inevitably evolves over time.  See [@winters2020software] for some excellent discussions of the principles of software reliability.
- Running biostatistical software on diseases with known biology and checking that this known biology is re-capitulated.  This check evaluates both: a) the overall reliability of the software and b) that I have correctly understood how to use it.

