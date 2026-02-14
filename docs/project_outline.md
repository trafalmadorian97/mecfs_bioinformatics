

# Project Outline
## Overview

The ME/CFS Biostatistics repo is a home for reproducible, open-source analysis of datasets pertaining to ME/CFS and other poorly-understood diseases.

## Principles

 There are a few key principles that govern how software should be written in this repo.



**Automated Reproducibility**:  Reproducibility is a central tenet of science, but a great deal of scientific data analysis is difficult to reproduce.  Even when the requisite data is public, reproducing published analysis on this data can be laborious.  In this repo, we mitigate this problem by automating reproducibility: **it should be possible to reproduce any of the main analyses with a few lines of Python code**. 



**Automated Testing**: Software engineering is ["programming integrated over time"](https://abseil.io/resources/swe-book/html/ch01.html#:~:text=Within%20Google%2C%20we%20sometimes%20say,software%20in%20the%20first%20place.)[@winters2020software]. Software engineering thus involves challenges unseen in short-lived programming projects.

One such challenge is that in any actively developed repository, the introduction of new features can break previously-developed features. A powerful mitigation against this risk is the use of automated unit and integration tests.  Automated tests can verify that old features continue to function as the codebase changes.  In this repo, we aim to **protect key features with automated testing.**  See[ Winters et al.](https://abseil.io/resources/swe-book/html/ch12.html)[@winters2020software] for discussion of the principles of software testing.



