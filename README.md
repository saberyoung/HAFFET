# SNDA(SN Data Analysis) to analyze SN photometric and spectroscopic data

This is the version of Docker set out to design for Windows users  by Marc Shen.

Install and Use
==

### Install the Docker

Get and install Docker from https://www.docker.com/get-started/

### Install the program

First, download this program.

Then do as follows.

```bash
mv -f arm.Dockerfile Dockerfile
# If the computer is arm(arm32, arm64 .etc) architecture.

mv -f amd.Dockerfile Dockerfile
# If the computer is amd(x86, x64 .etc) architecture.
```

```bash
sudo docker build -t test:v1 .
```

### Use

Run the program with gui.

```bash
docker run -ti --rm -e DISPLAY=host.docker.internal:0.0 test:v1 sdapy_gui
```

## Description of different systems

### Windows

#### WSL

You need WSL to run Docker.

#### XServer

To use this program you need to use this software together.

https://sourceforge.net/projects/vcxsrv/

### Linux & Mac

Under Linux, you need to adjust the `DISPLAY` for your computer system.

## Details

Linux and Mac users are advised to use this program with the master branch page.

**This method is mainly designed for Windows system.** 

Of course, Linux and Mac systems can also use this program, but it is not necessary.

## Q&A

This is a brand new test, you can give us feedback on specific questions.

