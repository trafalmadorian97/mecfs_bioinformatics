

# Project Outline
## Overview

The ME/CFS Biostatistics repo is a home for reproducible, open-source analysis of biostatistical datasets pertaining to ME/CFS and other poorly-understood diseases.

## Principles

 There are a few key principles that govern how software should be written in this repo.



**Automated Reproducibility**:  Reproducibility is a central tenet of science, but a great deal of modern scientific data analysis is difficult to reproduce.  Even when the requisite data is publicly available, reproducing published analysis on this data is often extremely labor-intensive.  In this repo, we seek to mitigate this problem by automating reproducibility: **it should be possible to reproduce any of the main analyses with a few lines of Python code**. 



**Automated Testing**: Titus Winters defines software engineering as ["programming integrated over time"](https://abseil.io/resources/swe-book/html/ch01.html#:~:text=Within%20Google%2C%20we%20sometimes%20say,software%20in%20the%20first%20place.). Software engineering thus involves challenges unseen in short-lived programming projects.

A central time-related software engineering risk is that in any actively developed code repository, the introduction of new features can break previously-developed features. One powerful mitigation against this risk is the use of automated unit and integration tests.  Automated tests can verify that old features continue to function as the codebase changes..  In this repo, we aim to **protect key features with automated testing.**  See[ Winters et al.](https://abseil.io/resources/swe-book/html/ch12.html)[@winters2020software] for discussion of the principles of software testing.



