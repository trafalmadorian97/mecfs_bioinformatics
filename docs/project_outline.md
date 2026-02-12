

# Project Outline
## Overview

The ME/CFS Biostatistics repo is a home for reproducible, open-source analysis of biostatistical datasets pertaining to ME/CFS and other poorly-understood diseases.

## Principles

 There are a few key principles that govern how software should be written in this repo.



**Automated Reproducibility**:  Reproducibility is a central tenet of science, but a great deal of modern scientific data analysis is difficult to reproduce.  Even when the requisite data is publicly available, reproducing published analysis on this data is often extremely labor-intensive.  In this repo, we seek to mitigate this problem by automating reproducibility: **it should be possible to reproduce any of the main analyses with a few lines of Python code**. 



**Automated Testing**: Titus Winters defines software engineering as ["programming integrated over time"](https://abseil.io/resources/swe-book/html/ch01.html#:~:text=Within%20Google%2C%20we%20sometimes%20say,software%20in%20the%20first%20place.). The discipline of software engineering thus involves a host of challenges that are not seen in short-lived programming projects.

A central time-related software engineering challenge is that in any actively developed code repository, the introduction of new functionality always risks breaking previous functionality. One powerful mitigation against this challenge is the use of unit and integration tests that verify that old functionality continues to work.  In this repo, we aim to **protect all key functionality with automated testing.**  See[ Winters et al.](https://abseil.io/resources/swe-book/html/ch12.html)[@winters2020software] for discussion of the principles of software testing.



[//]: # ()
[//]: # (**Verification**: As much as is feasible, code and analysis should be verified.)

[//]: # ()
[//]: # (Verification techniques include:)

[//]: # ()
[//]: # (- Running biostatistical software on diseases with known biology and checking that this known biology is re-capitulated.  This check evaluates both: a&#41; the reliability of the scientific assumptions underlying the software and b&#41; that we have correctly understood how to use the software.)

[//]: # ()
