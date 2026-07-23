

# Project Outline
## Overview

OnThe ME/CFS Biostatistics repo is a home for reproducible, open-source analysis of datasets pertaining to ME/CFS and other poorly-understood diseases.

## Principles

 There are a few key principles that govern how software should be written in this repo.



### Automated Reproducibility

Reproducibility is central to science, but a great deal of scientific data analysis is difficult to reproduce.  Even when the requisite data is public, reproducing published analysis can be laborious.  In this repo, we mitigate this problem by automating reproducibility: **it should be possible to reproduce any of the main analyses with a few lines of Python code**. 



### Automated Testing

Software engineering is [programming integrated over time](https://abseil.io/resources/swe-book/html/ch01.html#:~:text=Within%20Google%2C%20we%20sometimes%20say,software%20in%20the%20first%20place.)[@winters2020software]. It thus involves challenges unseen in short-lived programming projects.

One such challenge is the tendency of new features to break previously-developed features. A powerful mitigation against this risk is the use of automated unit and integration tests.  Automated tests can verify that old features continue to function as the codebase changes.  In this repo, we aim to **protect key features with automated testing.**  See[ Winters et al.](https://abseil.io/resources/swe-book/html/ch12.html)[@winters2020software] for discussion of the principles of software testing.



## AI Usage Rules

- **AI for coding**:  There can be no doubt that AI coding agents are productivity-enhancing.  However, current agents demonstrate what has been called "jagged intelligence"\footnote{The phrase was coined by Andrej Karpathy.}: they can exhibit superhuman capabilities on one task and then make obvious mistakes on another.  For this reason, **AI-written code must be reviewed and understood in detail by a human.** 
- **AI should not be used to generate write-ups for the `docs/` directory**: readers need to be able to trust that the discussions and analyses on the ME/CFS Bioinformatic site correspond to the thoughts of a human.