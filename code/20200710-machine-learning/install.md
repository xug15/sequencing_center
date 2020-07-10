# install.md
## Docker 安装
```sh
# Download tensorflow
docker pull tensorflow/tensorflow

# Start a CPU-only container
$ docker run -it --rm tensorflow/tensorflow bash

# Start a GPU container, using the Python interpreter.
$ docker run -it --rm --runtime=nvidia tensorflow/tensorflow:latest-gpu python

# Run a Jupyter notebook server with your own notebook directory (assumed here to be ~/notebooks). To use it, navigate to localhost:8888 in your browser.
$ docker run -it --rm -v $(realpath ~/notebooks):/tf/notebooks -p 8888:8888 tensorflow/tensorflow:latest-jupyter
# xugang computer example.
docker run -it --rm -v /Users/xugang/Documents/c-pycharm/git/sequencing_center/code/20200710-machine-learning/notebooks:/tf/notebooks -p 8888:8888 tensorflow/tensorflow:latest-jupyter


```
