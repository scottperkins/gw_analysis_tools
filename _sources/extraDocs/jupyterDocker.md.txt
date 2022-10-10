
Hosting jupyter servers in Docker containers can be extremely useful, ~~as this typically doesn't involve logging into the container at all,~~[^1] and you can just connect to your server in a browser.

[^1]: Can't figure this out right now

# Dockerfiles

A typical docker file for this purpose looks like

```docker
FROM debian:latest # or whatever base image you want

RUN pip3 install --upgrade pip
RUN pip3 install jupyter jupyterlab #plus whatever else

CMD ["jupyter-lab", "--no-browser", "--allow-root", "--ip=0.0.0.0", "--port=8989"] # The port doesn't have to be the same, just some random port not in use
```

# Running container

Now, with image in hand, navigate to the directory of interest and run the following command:

```bash
docker run -d --rm -p 8989:8989 -v $(PWD):/opt/project scottperkins/gwat.jupyter
```

where the port matches the port used in the image (at least the second one) and the "/opt/project" directory is whatever you want to call it (directory inside the container). The image (scottperkins/gwat.jupyter) is obviously the image of the container you want to run

# Logging in

For now, you have to log into the container with

```bash
docker exec -it <container-id> bash
```
 Then, inside the container, run

```bash
jupyter server list
```

Go to localhost:8989 in a browser and copy the token from the docker container to log in.

You're in! Now you can exit out of the interactive session of the container if you want.
