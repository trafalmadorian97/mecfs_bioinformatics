# Docker

Certain bioinformatic software packages are distributed as Docker images.  An example is [MiXeR](../Bioinformatics_Concepts/Mixer.md). Thus, running tasks that use these software packages ({{ api_link("MixerTask", "mecfs_bio.build_system.task.mixer.mixer_task") }} for example) requires installation of Docker.

- Follow the directions [here](https://docs.docker.com/engine/install/ubuntu/) to install Docker on Ubuntu Linux. There are similar instructions for installing on other operating systems
- To be able to run docker images without use of `sudo`, you need to be added to the Docker usergroup.  You can do this via

```
sudo usermod -aG docker $USER
```




