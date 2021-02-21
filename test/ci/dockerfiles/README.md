# HowTo with docker

## update an image and push it to the registry

1. go to folder with dockerfile, *e.g.* `cd test/ci/dockerfiles/ubuntu_18_04/`
2. build new image: `docker build -t registry.gitlab.com/uguca/uguca/ubuntu:18.04 .`
3. push to gitlab: `docker push registry.gitlab.com/uguca/uguca/ubuntu:18.04`

For more info, check: gitlab.com > packages & registries > container registry > CLI Commands


## test locally on your computer inside an image from registry

1. go to base folder of *uguca* `cd uguca`
2. `mkdir build-docker`
3. `docker run -it -u 1000:1000 -v $PWD:/uguca -w /uguca/build-docker registry.gitlab.com/uguca/uguca/ubuntu:18.04 bash` 

(`1000:1000` is assuming that you are first user of local computer)


## Set-up of docker software on computer

1. Install docker on linux computer
   * follow [https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository](https://docs.docker.com/engine/install/ubuntu/#install-using-the-repository)
   * `sudo apt-get remove docker docker-engine docker.io containerd runc`
   * `sudo apt-get update`
   * `sudo apt-get install apt-transport-https ca-certificates curl gnupg-agent software-properties-common`
   * `curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -`
   * `sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable` (for linux mint had to use `bionic` instead of `$(lsb_release -cs` check registry)
   * `sudo apt-get update`
   * `sudo apt-get install docker-ce docker-ce-cli containerd.io`
   * add yourself to docker group: `sudo gpasswd -a username docker` (logout and login)
2. Register to gitlab.com: `docker login registry.gitlab.com`
