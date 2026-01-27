
# Synapse

Some datasets must be downloaded from [Synapse](https://www.synapse.org/), and require a Synapse access token.

## Steps to Set Up Synapse Access Token


1.  Go to https://www.synapse.org/, and register to create an account.
2. Click on Account Settings -> Personal Access Tokens -> Manage Personal Access Tokens -> Create New Token
3. Create a file `~/.synapseConfig` at the root of your user directory.  Edit it so it has the following contents

```bash
[default]
authtoken = <YOUR_AUTH_TOKEN>
```

