
# Synapse

Some datasets must be downloaded from [Synapse](https://www.synapse.org/). An example is the [UK Biobank Pharma Proteomics Project Summary Statistics](../Data_Sources/UKBB_PPP.md)  While the build system can automate these download, each user must provide their own Synapse Personal Access Token.

## Steps to Set Up Synapse Access Token


1.  Go to https://www.synapse.org/, and register to create an account.
2. Click on Account Settings -> Personal Access Tokens -> Manage Personal Access Tokens -> Create New Token
3. Create a file `~/.synapseConfig` at the root of your user directory.  Edit it so it has the following contents

```bash
[default]
authtoken = <YOUR_AUTH_TOKEN>
```

