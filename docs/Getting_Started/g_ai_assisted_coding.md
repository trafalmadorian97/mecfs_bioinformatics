# AI Usage

As of 2026, AI-assisted coding is new and rapidly-changing, so there are not a canonical body of best-practices.  Below I:


- Attempt to formulate a set of rules for AI usage that enable us to take advantage of its AI's productivity benefits without sacrificing correctness of misleading readers; 
- Describe some practical advice I have gathered from my limited experience.

## Rules

- **AI for coding**:  There can be no doubt that AI coding agents enhance productivity.  However, current agents demonstrate what has been called "jagged intelligence"[^foot]: they can exhibit superhuman capabilities on one task and then make obvious mistakes on another.  For this reason, **AI-written code must be reviewed and understood in by a human.**
- **AI should not be used to generate write-ups for the `docs/` directory**.  Readers need to be able to trust that the discussions and analyses on the ME/CFS Bioinformatics site correspond to the considered thoughts of a human.



## Advice
- AI agents work best when they have some verifiable source-of-truth against which they can check their work.  It is often worthwhile to start an AI coding session by first getting the agent to create an automated test for verification, then asking the agent to write its code to make the test pass.  This is consistent with the classical software engineering practice of [test-driven-development](https://en.wikipedia.org/wiki/Test-driven_development)
- I have noticed that AI agents have a tendency to ignore existing repo infrastructure and re-invent the wheel.  For instance, they may fail to find an already-existing utility function that accomplishes a task, and re-implement something very similar. To avoid this, it can be helpful to point out relevant infrastructure before initiating  a task.
- AI agents tend to write bad comments.  These comments can be bad because they are written in hard-to-decipher "LLM-speak", or because they spend many paragraphs elaborating on an obscure point that is irrelevant to the general reader.  The only solution I have found is to manually edit comments.



[^foot]: The phrase was coined by Andrej Karpathy.
